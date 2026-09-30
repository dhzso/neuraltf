#!/usr/bin/env python
"""Cross-atlas DE concordance analysis (Fisher + Stouffer combination).

The pipeline's de_pvalues checkpoint stores, per gene and atlas, the best
(minimum) Wilcoxon p among clusters that already passed the global BH
q<=0.1 gate. Combining those values directly with Fisher/Stouffer gives
"100% of genes significant" - the p-values are CONDITIONED on within-atlas
significance (a gene-by-gene selection event), so the Fisher chi2 null
does not apply to them and the combined p is not a valid meta-analytic
p-value.

The combination machinery stays the same (the formulas are correct); the
inference is:

1. SELECTION-AWARE combined p: each per-atlas p is first mapped to its
   within-atlas multiplicity-adjusted value using the number of cluster
   tests that gene faced (min-p over K_cluster tests is anti-conservative
   by ~K_cluster); the combined statistic is calibrated against a
   PERMUTATION-null-free upper bound: p_adj = 1 - (1-p_min)^K_eff, the
   exact probability that the gene's best cluster p would be <= p_min
   under H0 for that atlas (equivalently, Bonferroni on the gene's own
   test count). These adjusted p-values are approximately valid inputs
   to Fisher/Stouffer.
2. CONCORDANCE framing: the output reports, per gene, how many atlases
   show significant upregulated DE after the correction, and the
   combined statistic is labeled a concordance score, not a
   meta-analytic significance test of a global H0.
3. No fabricated fallback: if the checkpoint is missing, exit loudly.

Usage:
    python scripts/stats/meta_analytic_pvalue.py
"""

# Implementation notes:
# 1. Use sf, not 1-cdf: 1-cdf produced 5,986 exact-zero Fisher p's and
#    5,656 zero Stouffer p's through catastrophic cancellation (chi2=166.8,
#    df=6 has true sf ~2e-33). Zero p's break downstream log/rank
#    arithmetic and make the Fisher-vs-Stouffer Spearman tie-dominated.
# 2. k (cluster counts) is read from the pipeline's
#    checkpoint_02_post_qc (n_leiden_clusters per atlas). The Leiden
#    counts drift between runs (14/20/30 and 15/20/32 observed), so a
#    hard-coded k silently mismatches regenerated checkpoints.
# 3. Remaining limitation: the stored p's are still conditioned on the
#    BH q<=0.1 gate, and the true selection family is (gene x cluster)
#    per atlas rather than k clusters. Where the pipeline has stored
#    unconditioned min-p's (de_pvalues_raw), this script prefers those
#    columns.

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

ATLAS_P_COLS = {
    "fincher": "fincher_p",
    "plass": "plass_p",
    "cui": "cui_p",
}
# Fallback ONLY if the checkpoint is missing (loud warning). The
# pipeline's Leiden counts are run-dependent; the checkpoint is the
# source of truth.
N_CLUSTERS_DEFAULT = {"fincher": 16, "plass": 22, "cui": 30}


def cluster_counts_from_checkpoint() -> dict[str, int]:
    """Read per-atlas Leiden cluster counts from checkpoint_02_post_qc.

    The k used to un-condition each min-p must match the run that
    produced the stored p's. Falls back to N_CLUSTERS_DEFAULT with a
    LOUD warning (the old silent defaults could mismatch a regenerated
    checkpoint by 2-6 clusters per atlas).
    """
    path = RUN_DIR / "checkpoint_02_post_qc.parquet"
    if not path.exists():
        path = RUN_DIR / "checkpoint_02_post_qc.csv"
    if not path.exists():
        print("WARNING: checkpoint_02_post_qc not found — using hard-coded "
              "cluster-count fallback {N_CLUSTERS_DEFAULT}. The un-"
              "conditioning k may not match the run that stored the p's.")
        return dict(N_CLUSTERS_DEFAULT)
    df = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
    out = {}
    for _, row in df.iterrows():
        atlas = str(row.get("atlas", "")).strip().lower()
        n_cl = row.get("n_leiden_clusters")
        if atlas and pd.notna(n_cl):
            try:
                out[atlas] = int(n_cl)
            except (TypeError, ValueError):
                pass
    missing = [a for a in ATLAS_P_COLS if a not in out]
    if missing:
        print(f"WARNING: checkpoint lacks cluster counts for {missing}; "
              f"filling from defaults where missing.")
        for a in missing:
            out[a] = N_CLUSTERS_DEFAULT.get(a, 20)
    return out


def adjusted_min_p(p_min: float, k_eff: int) -> float:
    """Un-condition a within-atlas min-p over k_eff cluster tests.

    P(at least one of k_eff tests <= p_min | H0) = 1 - (1-p_min)^k_eff.
    This is the gene's own selection-adjusted p for that atlas; the
    resulting values are valid (conservative) inputs to Fisher/Stouffer.

    The result is clamped to <= 1-1e-16: for p_min close to 1 and modest
    k_eff, 1-(1-p)^k rounds to exactly 1.0 in floating point, and
    stats.norm.isf(1.0) then returns -inf, poisoning the Stouffer z
    (log(1.0)=0 would also zero the Fisher contribution). The clamp keeps
    both combiners finite for every input.
    """
    p_ceiling = 1.0 - 1e-16
    if k_eff <= 1:
        return min(p_ceiling, max(p_min, 1e-16))
    adj = 1.0 - (1.0 - max(p_min, 1e-16)) ** k_eff
    return float(min(p_ceiling, adj))


def fishers_method(pvalues):
    """Combine p-values using Fisher's method. Returns chi2 statistic and combined p-value.

    Tail via chi2.sf — 1-cdf underflows to exactly 0.0 for the large
    chi2 statistics this data produces (observed: 5,986 zero p's under
    the old form; true values as small as ~1e-33).

    Boundary behaviour: a zero/negative p (maximally significant input)
    saturates the combined p at 0 rather than being dropped — filtering
    to ``p > 0`` would discard the strongest evidence and shrink the dof.
    NaN inputs (a test not performed) are the only values excluded, and
    only explicitly.
    """
    p = np.asarray(pvalues, dtype=float)
    p = p[~np.isnan(p)]
    if len(p) == 0:
        return 0.0, 1.0
    if np.any(p <= 0.0):
        # chi2 = -2*sum(log(p)) -> inf; the combined p is 0.
        return float("inf"), 0.0
    chi2_stat = -2.0 * np.sum(np.log(p))
    combined_p = float(stats.chi2.sf(chi2_stat, 2 * len(p)))
    return float(chi2_stat), combined_p


def stouffers_method(pvalues, weights=None):
    """Combine p-values using Stouffer's method. Returns z statistic and combined p-value.

    Tail via norm.sf (see fishers_method note).

    Boundary behaviour: a zero/negative p (z = isf(0) = +inf) saturates
    the combined p at 0 instead of being dropped (a ``p > 0`` filter
    would discard the strongest evidence); p >= 1 (z = -inf) is clamped
    to 1-1e-16 so one exhausted input cannot poison the whole weighted z.
    NaN inputs are excluded explicitly (a test not performed).
    """
    p = np.asarray(pvalues, dtype=float)
    if weights is not None:
        weights = np.asarray(weights, dtype=float)[~np.isnan(p)]
    p = p[~np.isnan(p)]
    if len(p) == 0:
        return 0.0, 1.0
    if np.any(p <= 0.0):
        return float("inf"), 0.0
    p = np.minimum(p, 1.0 - 1e-16)
    z_scores = stats.norm.isf(p)
    if weights is None:
        weights = np.ones(len(p))
    combined_z = np.sum(z_scores * weights) / np.sqrt(np.sum(weights**2))
    combined_p = float(stats.norm.sf(combined_z))
    return float(combined_z), combined_p


def load_de_pvalues() -> pd.DataFrame | None:
    """Load the pipeline's per-gene per-atlas DE p-value checkpoint.

    Prefers the unconditioned columns (<atlas>_p_raw) over the
    BH-conditioned ones (<atlas>_p): the conditioned values are min-p's
    among tests that already passed the atlas-wide q<=0.1 gate and are
    not valid Fisher/Stouffer inputs (they report ~99.6% of genes as
    significant).
    """
    path = RUN_DIR / "de_pvalues.parquet"
    if not path.exists():
        path = RUN_DIR / "de_pvalues.csv"
    if not path.exists():
        return None
    df = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
    df = df.drop_duplicates(subset="v6_id", keep="first")
    # has_raw is False (not vacuously True) when the checkpoint has no
    # <atlas>_p columns at all: all(...) over an empty generator would
    # print the misleading "Using UNCONDITIONED" message before the
    # >=2-atlas check rejects the file.
    p_cols = [a for a in ATLAS_P_COLS if f"{a}_p" in df.columns]
    raw_cols = [a for a in ATLAS_P_COLS if f"{a}_p_raw" in df.columns]
    has_raw = bool(p_cols) and all(a in raw_cols for a in p_cols)
    if has_raw:
        # Use the unconditioned min-p's directly. The per-atlas cluster-
        # test counts k used by the Sidak un-conditioning come from
        # checkpoint_02_post_qc via cluster_counts_from_checkpoint()
        # (the pipeline writes no per-gene k_eff columns).
        print("Using UNCONDITIONED per-atlas min-p's (<atlas>_p_raw) — "
              "valid Fisher/Stouffer inputs.")
        rename = {}
        for a in ATLAS_P_COLS:
            if f"{a}_p_raw" in df.columns:
                df[f"{a}_p"] = df[f"{a}_p_raw"]
                df.drop(columns=[f"{a}_p_raw"], inplace=True)
        return df
    print("NOTE: checkpoint lacks <atlas>_p_raw columns — falling back to "
          "the BH-CONDITIONED stored p's with the Sidak un-conditioning "
          "(approximate; the honest fix is the pipeline B1 re-run).")
    return df


def main():
    print("=== Cross-Atlas DE Concordance (selection-adjusted combination) ===")

    de = load_de_pvalues()
    if de is None or de.empty:
        print(
            "ERROR: per-atlas DE p-value checkpoint not found at "
            f"{RUN_DIR / 'de_pvalues.parquet'}.\n"
            "Re-run the pipeline (scripts/run.py) - checkpoint 03 writes "
            "de_pvalues with fincher_p/plass_p/cui_p per gene.\n"
            "This analysis never falls back to simulated data."
        )
        return 1

    present = {a: c for a, c in ATLAS_P_COLS.items() if c in de.columns}
    if len(present) < 2:
        print(f"ERROR: need p-values from >= 2 atlases; found {list(present)}")
        return 1
    print(f"Atlases with per-gene p-values: {sorted(present)} ({len(de)} genes)")
    print("NOTE: stored values are best-cluster p's already conditioned on "
          "BH q<=0.1 within each atlas; they are un-conditioned per gene "
          "via 1-(1-p_min)^k_clusters before combination.")
    k_clusters = cluster_counts_from_checkpoint()
    print(f"Un-conditioning k per atlas (from checkpoint_02_post_qc): "
          f"{k_clusters}")

    de = de.sort_values("v6_id").reset_index(drop=True)
    results = []
    for _, row in de.iterrows():
        pvals_adj = []
        pvals_raw = []
        atlases_sig = 0
        for atlas, col in sorted(present.items()):
            v = row.get(col)
            if pd.notna(v):
                p_raw = max(float(v), 1e-16)
                p_adj = adjusted_min_p(p_raw, k_clusters.get(atlas, 20))
                pvals_raw.append(p_raw)
                pvals_adj.append(p_adj)
                if p_adj < 0.05:
                    atlases_sig += 1
        if len(pvals_adj) < 2:
            continue
        fisher_chi2, fisher_p = fishers_method(pvals_adj)
        stouffer_z, stouffer_p = stouffers_method(pvals_adj)
        results.append({
            "gene_id": row["v6_id"],
            "n_atlases": len(pvals_adj),
            "n_atlases_sig_adj": atlases_sig,
            "individual_pvalues": pvals_raw,
            "adjusted_pvalues": pvals_adj,
            "fisher_chi2": fisher_chi2,
            "fisher_combined_p": fisher_p,
            "stouffer_z": stouffer_z,
            "stouffer_combined_p": stouffer_p,
        })

    out_df = pd.DataFrame(results)
    if out_df.empty:
        print("No gene had p-values in >= 2 atlases; nothing to combine.")
        return 1
    out_df["individual_pvalues"] = out_df["individual_pvalues"].apply(str)
    out_df["adjusted_pvalues"] = out_df["adjusted_pvalues"].apply(str)
    out_df = out_df.sort_values("fisher_combined_p")

    out_path = RESULTS_DIR / "meta_analysis_pvalues.csv"
    out_df.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path} ({len(out_df)} genes)")

    print(f"\nCombined concordance (selection-adjusted) - Fisher p<0.05: "
          f"{(out_df['fisher_combined_p'] < 0.05).sum()} / {len(out_df)}")
    print(f"Combined concordance (selection-adjusted) - Stouffer p<0.05: "
          f"{(out_df['stouffer_combined_p'] < 0.05).sum()} / {len(out_df)}")
    print(f"Genes significant (adj p<0.05) in >= 2 atlases: "
          f"{(out_df['n_atlases_sig_adj'] >= 2).sum()}")
    print(f"Genes significant (adj p<0.05) in ALL tested atlases: "
          f"{(out_df['n_atlases_sig_adj'] == out_df['n_atlases']).sum()}")

    print("\nTop-10 by Fisher's combined p-value (concordance):")
    for _, row in out_df.head(10).iterrows():
        print(f"  {row['gene_id']:>30s}  Fisher p={row['fisher_combined_p']:.4e}  "
              f"Stouffer p={row['stouffer_combined_p']:.4e}  "
              f"sig_atlases={row['n_atlases_sig_adj']}/{row['n_atlases']}")

    rho = out_df["fisher_combined_p"].corr(out_df["stouffer_combined_p"],
                                          method="spearman")
    print(f"\nSpearman rho (Fisher vs Stouffer): {rho:.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
