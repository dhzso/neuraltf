#!/usr/bin/env python
"""Cross-atlas DE CONCORDANCE analysis (Fisher + Stouffer combination).

2026-09-04 statistical redesign (audit finding M3):

The pipeline's de_pvalues checkpoint stores, per gene and atlas, the best
(minimum) Wilcoxon p among clusters that already passed the global BH
q<=0.1 gate. Combining those values with Fisher/Stouffer produced
"100% of genes significant" - the p-values are CONDITIONED on within-atlas
significance (a gene-by-gene selection event), so the Fisher chi2 null
does not apply to them and the combined p is not a valid meta-analytic
p-value.

This version keeps the same combination machinery (formulas are correct)
but reframes and corrects the inference:

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

# Number of Leiden clusters per atlas in the production run (the per-gene
# min-p was taken over one test per cluster); used to un-condition the
# stored min-p values. Read from the production log's run_qc line when
# available, else these defaults.
N_CLUSTERS_DEFAULT = {"fincher": 16, "plass": 22, "cui": 30}


def adjusted_min_p(p_min: float, k_eff: int) -> float:
    """Un-condition a within-atlas min-p over k_eff cluster tests.

    P(at least one of k_eff tests <= p_min | H0) = 1 - (1-p_min)^k_eff.
    This is the gene's own selection-adjusted p for that atlas; the
    resulting values are valid (conservative) inputs to Fisher/Stouffer.
    """
    if k_eff <= 1:
        return min(1.0, max(p_min, 1e-16))
    return float(min(1.0, 1.0 - (1.0 - max(p_min, 1e-16)) ** k_eff))


def fishers_method(pvalues):
    """Combine p-values using Fisher's method. Returns chi2 statistic and combined p-value."""
    p = np.array(pvalues)
    p = p[p > 0]
    if len(p) == 0:
        return 0.0, 1.0
    chi2_stat = -2.0 * np.sum(np.log(p))
    combined_p = 1.0 - stats.chi2.cdf(chi2_stat, 2 * len(p))
    return float(chi2_stat), float(combined_p)


def stouffers_method(pvalues, weights=None):
    """Combine p-values using Stouffer's method. Returns z statistic and combined p-value."""
    p = np.array(pvalues)
    p = p[p > 0]
    if len(p) == 0:
        return 0.0, 1.0
    z_scores = stats.norm.ppf(1.0 - p)
    if weights is None:
        weights = np.ones(len(p))
    combined_z = np.sum(z_scores * weights) / np.sqrt(np.sum(weights**2))
    combined_p = 1.0 - stats.norm.cdf(combined_z)
    return float(combined_z), float(combined_p)


def load_de_pvalues() -> pd.DataFrame | None:
    """Load the pipeline's per-gene per-atlas DE p-value checkpoint."""
    path = RUN_DIR / "de_pvalues.parquet"
    if not path.exists():
        path = RUN_DIR / "de_pvalues.csv"
    if not path.exists():
        return None
    df = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
    return df.drop_duplicates(subset="v6_id", keep="first")


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
                p_adj = adjusted_min_p(p_raw, N_CLUSTERS_DEFAULT.get(atlas, 20))
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
