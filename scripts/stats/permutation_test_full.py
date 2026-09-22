#!/usr/bin/env python
"""Permutation test on the SAME multi-stream statistic the pipeline ranks.

2026-09-04 FULL REDESIGN:

The previous null was degenerate: it permuted only Fincher/Plass labels,
silently RETAINED the observed (King/Cui-inflated) expression for genes
with no permuted DE hit, forced null specificity to 1.0 via a dead regex,
and ran on 2,000-cell TF-only re-clustered atlases. 99.6% of genes landed
at the floor p-value - the test measured nothing.

Exchangeable-null design (this version):
- Leiden clustering is run ONCE per atlas in ``prepare_atlas()`` to
  establish the same cluster structure the real scores were computed on.
  Each permutation draw SHUFFLES the resulting cluster labels (a random
  partition of equal cluster sizes), then re-runs Wilcoxon DE under the
  permuted assignment — Leiden itself is NOT re-run per draw.
- Only 2 of 11 evidence streams — EXPRESSION and SPECIFICITY — are
  recomputed from the permuted atlas DE.  The remaining 9 streams
  (reproducibility, rnai, correlation, neural_enriched,
  neural_specificity, perez_lineage, perez_influence, fincher_brain,
  cui_temporal) are HELD FIXED at their observed values in both the real
  and null statistics because permuting cluster labels cannot change an
  external table or a non-DE-derived quantity.
- The King expression floor min(1, max_log2FC/5) and King fractional
  breadth specificity floor 1-(n-1)/(N-1) enter BOTH real and null via
  max(atlas value, King floor), exactly as the pipeline integrates them.
  Consequently a gene whose score is fully explained by table evidence
  (e.g. King log2FC >= 5 saturating expression at 1.0) has null == real
  by construction: its p-value is set to 1.0 and it is flagged
  untestable_by_permutation=True — the honest statement that its rank
  rests on external tables, not on cluster-specific atlas expression.
- Stream PRESENCE mirrors the pipeline: expression/specificity are
  present in the null iff the gene has permuted atlas evidence OR a King
  floor; otherwise the streams are absent (NaN) as in the real score.
- Empirical p = (#{null >= real} + 1) / (n + 1), one-tailed add-one,
  computed ONLY for testable genes (those whose real score exceeds the
  label-independent ceiling by more than floating-point tolerance).
  Untestable genes receive p = 1.0 by construction.
- BH-FDR correction is applied ONLY over testable genes.  Untestable
  genes receive q = NaN (they were never tested).  This prevents
  degenerate floor p-values from inflating the tested family's q.

The null answers: "could this gene's integrated score arise from random
cluster structure in the atlases, given its label-independent evidence?"

LIMITATIONS AND EXCHANGEABILITY GAPS:
- Only expression and specificity (2/11 streams) are under the null.
  The p-value is a statement about the ATLAS-DE contribution only, not
  about the full integrated score.
- reproducibility (cross-atlas membership), fincher_brain, and
  cui_temporal are atlas-derived quantities but are held fixed.  A
  fully exchangeable null would permute these too, but that would
  require re-running the full pipeline per draw.
- The label-independent "ceiling" is NOT a strict upper bound: after
  renormalization over present streams, permuted denominators can
  shift slightly, producing surplus < 0 for some genes.  Genes whose
  surplus is within floating-point tolerance of zero are classified
  as untestable.
- In --candidates mode, BH is applied over the 143-gene neural panel
  (not the full 11,696-gene universe), which is more lenient than the
  pipeline's full-universe gate.

2026-09-06 --candidates mode (targeted resolution):
The production n=30 full-universe run is resolution-floored at p=1/31
(96.8% of genes pinned at the floor — no per-gene claim possible). The
--candidates flag restricts the test to the neural-candidate family
(rank_neural.csv, 143 genes) so n=1000 draws resolve p ~ 1e-3 with honest
BH WITHIN that family (Bonferroni alpha = 0.05/143 = 3.5e-4 needs
n>=2,860 for full family control; n=1000 resolves BH-q at the 143-gene
scale). The permutation machinery is unchanged — same atlas permutation,
same BH-FDR DE gate, same King floors and testability flags; only the
p-value aggregation family shrinks to the genes that matter for the
shortlist.

2026-09-21 AUDIT FIX (float-tie artifacts + BH family):
The previous version used >= counting for ALL genes, so untestable genes
whose real score was within floating-point epsilon of the ceiling received
the floor p-value (1/1001 ≈ 1e-3) instead of the correct p=1.0.  This
produced 12 false-significant genes (out of 74 marked significant),
including two in the top-10 shortlist (dd38342, dd16466).  BH was also
computed over all 143 genes, including untestable ones, deflating q for
everyone.  Fixed: untestable genes get p=1.0 by enforcement; BH runs
only over testable genes; untestable genes get q=NaN.

Usage:
    python scripts/stats/permutation_test_full.py --n-perm 1000 --candidates
    python scripts/stats/permutation_test_full.py --n-perm 30
"""

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
DATA_DIR = REPO / "projects" / "NeuralTF" / "data"
RAW_DIR = REPO / "datasets" / "raw"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

STREAMS = ["expression", "specificity", "reproducibility", "rnai",
           "correlation", "neural_enriched", "neural_specificity",
           "perez_lineage", "perez_influence", "fincher_brain", "cui_temporal"]
# 2026-09-19: import the single-source-of-truth weights (previously a
# hand-copied vector).
import os
sys.path.insert(0, os.environ.get("BIOFORGE_SRC", str(REPO / "src")))
from bioforge.evidence.scoring import DEFAULT_WEIGHTS as _DW  # noqa: E402
W_BY_NAME = {getattr(k, "value", k): float(v) for k, v in _DW.items()}
W_DEFAULT = np.array([W_BY_NAME[s] for s in STREAMS])
assert abs(W_DEFAULT.sum() - 1.0) < 1e-9

FINCHER_PATH = REPO / "datasets" / "processed" / "fincher_subsample.h5ad"
PLASS_PATH = REPO / "datasets" / "processed" / "plass_v6.h5ad"
CUI_PATH = REPO / "datasets" / "processed" / "cui_v6.h5ad"
BRIDGE_PATH = DATA_DIR / "bridge.csv"
KING_MMC4 = (RAW_DIR / "Supplementary_Data_ King_2024"
             / "1-s2.0-S2211124724001712-mmc4.xlsx")
KING_ATLAS_TSV = DATA_DIR / "king_atlas.tsv"
MASTER_CATALOG = DATA_DIR / "master_tf_catalog.csv"

FDR_THRESHOLD = 0.1
L2FC_EPS = 1e-9
EXPR_CAP = 5.0


def load_tf_seed() -> set[str]:
    """King mmc4 TFs expanded with the master catalog (Perez MOESM5),
    exactly like pipeline.load_reference_tables."""
    df = pd.read_excel(KING_MMC4, sheet_name="TF")
    tf_ids = set(df.loc[df["TF?"].notna(), "Gene ID"].astype(str))
    tf_ids |= {tid[:-2] for tid in tf_ids if tid.endswith("_1")}
    if MASTER_CATALOG.exists():
        master = pd.read_csv(MASTER_CATALOG, dtype=str)
        if "v6_id" in master.columns:
            tf_ids |= set(master["v6_id"].dropna().str.strip()) - {"", "nan"}
    return tf_ids


def load_tf_seed_king_only() -> set[str]:
    """King mmc4-only seed - the exact set the pipeline forces into HVGs
    in run_qc (self.tf_ids_king)."""
    df = pd.read_excel(KING_MMC4, sheet_name="TF")
    return set(df.loc[df["TF?"].notna(), "Gene ID"].astype(str))


def load_bridge():
    from bioforge.evidence import load_bridge as _load_bridge
    return _load_bridge(BRIDGE_PATH)


def load_king_floors() -> tuple[dict[str, float], dict[str, float]]:
    """Label-INDEPENDENT King mmc7 floors, identical to the pipeline's
    gene-level push in integrate_king_atlas:

    expr floor  = min(1, max log2FC over ALL subcluster hits / 5)
    spec floor   = 1 - (nsub - 1)/(N_total - 1)   [fractional breadth,
                   (compartment, subcluster) pairs on both sides]

    Returns (expr_floor, spec_floor) keyed by v6_id. These enter BOTH the
    real and null statistics (max with the atlas component), so a gene
    saturated by the King table is honestly untestable by permutation.
    """
    if not KING_ATLAS_TSV.exists():
        return {}, {}
    king = pd.read_csv(KING_ATLAS_TSV, sep="\t")
    n_total = max(king.groupby(["compartment", "subcluster"]).ngroups, 2)

    expr_floor: dict[str, float] = {}
    spec_floor: dict[str, float] = {}
    for v6, hits in king.groupby("v6_id"):
        fcm = float(hits["log2fc"].max())
        expr_floor[str(v6)] = min(1.0, fcm / EXPR_CAP)
        nsub = hits.groupby(["compartment", "subcluster"]).ngroups
        spec_floor[str(v6)] = max(0.0, min(1.0, 1.0 - (nsub - 1) / (n_total - 1)))
    return expr_floor, spec_floor


def prepare_atlas(adata, name: str, tf_ids_king: set[str]):
    """Replicate the pipeline's run_qc ONCE per atlas (filter, normalize,
    log1p, HVG+forced King TFs, PCA, neighbors, igraph leiden 0.5,
    random_state fixed) to establish the cluster structure the real scores
    were computed on.

    NOTE: Leiden clustering is run ONCE here.  The permutation loop in
    main() then SHUFFLES the resulting cluster labels (via
    ``rng.permutation(obs['leiden'])``); it does NOT re-run Leiden per
    draw.  The null model is therefore a random partition of equal cluster
    sizes, not a re-clustering null.
    """
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.filter_cells(adata, min_counts=1)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=5000, batch_key=None)
        tf_in = [v for v in adata.var_names if v in tf_ids_king]
        adata.var.loc[tf_in, "highly_variable"] = True
        adata.raw = adata
        hvg = adata[:, adata.var.highly_variable].copy()
        sc.pp.pca(hvg, n_comps=50)
        sc.pp.neighbors(hvg, n_neighbors=10, n_pcs=40)
        sc.tl.leiden(hvg, resolution=0.5, flavor="igraph", n_iterations=2,
                     directed=False, random_state=42)
        adata.obs["leiden"] = hvg.obs["leiden"].values
    return adata


def permuted_atlas_streams(adata, atlas_name: str, tf_ids: set[str],
                           bridge, rng) -> dict:
    """One permutation: permute leiden labels, rerun Wilcoxon DE on the
    TF gene panel, and return {v6_id: (expr, spec)} derived from the
    permuted structure alone.

    expr = min(1, true_log2FC_best_cluster / 5)   (one-tailed)
    spec = 1 / n_sig_clusters                     (permuted breadth)
    """
    a = adata
    perm_labels = rng.permutation(a.obs["leiden"].astype(str).values)
    a.obs["leiden_perm"] = pd.Categorical(perm_labels)

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(a, "leiden_perm", method="wilcoxon")
    result = a.uns["rank_genes_groups"]
    clusters = list(result["names"].dtype.names)

    # Global BH over genes x clusters (matches pipeline._score_one_atlas)
    from statsmodels.stats.multitest import multipletests
    all_p, keys, idxs = [], [], []
    for cl in clusters:
        names = [str(g) for g in result["names"][cl]]
        for i, (g, p) in enumerate(zip(names, result["pvals"][cl])):
            all_p.append(float(p))
            keys.append((g, cl))
            idxs.append(i)

    _, qvals, _, _ = multipletests(all_p, alpha=FDR_THRESHOLD, method="fdr_bh")

    gene_best: dict[str, tuple[float, str]] = {}
    gene_sig: dict[str, set] = {}
    for (g, cl), q, i in zip(keys, qvals, idxs):
        lfc = float(result["logfoldchanges"][cl][i])
        pos = max(0.0, lfc)
        if q <= FDR_THRESHOLD and pos > 0:
            gene_sig.setdefault(g, set()).add(cl)
            if g not in gene_best or pos > gene_best[g][0]:
                gene_best[g] = (pos, cl)

    scores: dict[str, tuple[float, float]] = {}
    for gene, (pos_lfc, cl) in gene_best.items():
        if atlas_name == "fincher":
            v6 = bridge.v4_to_v6(gene)
        else:
            v6 = gene if gene in tf_ids else (gene + "_1" if gene + "_1" in tf_ids else None)
        if not v6:
            continue
        true_l2fc = _true_cluster_log2fc(a, gene, cl)
        n_sig = len(gene_sig.get(gene, set()))
        scores[v6] = (min(1.0, max(0.0, true_l2fc) / EXPR_CAP),
                      1.0 / n_sig if n_sig > 0 else 0.0)
    return scores


def _true_cluster_log2fc(adata, gene: str, cluster: str) -> float:
    """True log2FC from linear-space cluster means (same as the pipeline's
    _cluster_log2fc, on the permuted labels)."""
    try:
        gene_idx = adata.var_names.get_loc(gene)
    except (KeyError, ValueError):
        return 0.0
    labels = adata.obs["leiden_perm"].astype(str).values
    src = adata.raw if adata.raw is not None else adata
    col = src.X[:, gene_idx]
    vals = np.asarray(col.todense()).ravel() if hasattr(col, "todense") \
        else np.asarray(col).ravel()
    lin = np.expm1(vals.astype(np.float64))
    lin = np.where(lin < 0, 0.0, lin)
    in_cl = labels == str(cluster)
    n_in, n_out = int(in_cl.sum()), int((~in_cl).sum())
    if n_in == 0 or n_out == 0:
        return 0.0
    fc = (lin[in_cl].mean() + L2FC_EPS) / (lin[~in_cl].mean() + L2FC_EPS)
    return float(np.log2(fc)) if fc > 0 else 0.0


def integrated_score_with_renorm(S: np.ndarray, W: np.ndarray) -> float:
    """Compute integrated score with missing-data renormalization
    (identical to EvidenceScorer)."""
    mask = ~np.isnan(S)
    if not mask.any():
        return 0.0
    num = np.where(np.isnan(S), 0.0, S) @ W
    den = W[mask].sum()
    return num / den if den > 0 else 0.0


def main():
    parser = argparse.ArgumentParser(
        description="Permutation test on the full multi-stream statistic"
    )
    parser.add_argument("--n-perm", type=int, default=30)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--candidates", action="store_true",
                        help="restrict p-value aggregation to the neural-"
                             "candidate family (rank_neural.csv) — the "
                             "targeted high-resolution mode (n=1000 "
                             "recommended); permutation machinery unchanged")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    print(f"=== Permutation Test (exchangeable 3-atlas null, n={args.n_perm}) ===")

    real_rank_path = RUN_DIR / "rank.csv"
    if not real_rank_path.exists():
        print("Real rank.csv not found; run pipeline first")
        return 1
    real_rank = pd.read_csv(real_rank_path).drop_duplicates(
        subset="gene_id", keep="first")
    real_scores = dict(zip(real_rank["gene_id"], real_rank["integrated_score"]))
    print(f"Real candidates: {len(real_rank)}")

    # Observed per-gene stream vectors: label-INDEPENDENT streams stay
    # fixed under the null; expression/specificity are recomputed from
    # permuted atlases + the King floors.
    stream_cols = [c for c in STREAMS if c in real_rank.columns]
    stream_idx = {c: STREAMS.index(c) for c in stream_cols}
    observed = {}
    for _, row in real_rank.iterrows():
        S = np.full(len(STREAMS), np.nan)
        for c in stream_cols:
            v = row[c]
            S[stream_idx[c]] = float(v) if pd.notna(v) else np.nan
        observed[row["gene_id"]] = S

    # Label-independent King floors (enter BOTH real and null via max).
    king_expr_floor, king_spec_floor = load_king_floors()
    print(f"King floors: {len(king_expr_floor)} genes "
          f"({sum(1 for v in king_expr_floor.values() if v >= 1.0)} saturated at 1.0)")

    # Testability: a gene whose real score is fully explained by
    # label-independent evidence (King floors on expression/specificity)
    # cannot be distinguished from chance by cluster permutation - the
    # null's max(atlas_perm, King floor) >= King floor alone >= its real
    # contribution. The ceiling strips ALL atlas-derived expression and
    # specificity, keeping the King floor only where it exists.
    label_independent = {}
    for gid, S_obs in observed.items():
        S = S_obs.copy()
        ke = king_expr_floor.get(gid)
        ks = king_spec_floor.get(gid)
        S[0] = ke if ke is not None else np.nan   # strip atlas part
        S[1] = ks if ks is not None else np.nan
        label_independent[gid] = integrated_score_with_renorm(S, W_DEFAULT)
    n_untestable = sum(
        1 for gid in real_scores
        if real_scores[gid] <= label_independent.get(gid, 0.0) + 1e-9
    )
    print(f"Genes fully explained by label-independent evidence "
          f"(untestable by permutation): {n_untestable}")

    print("Loading atlases at production scale...")
    tf_ids = load_tf_seed()
    tf_ids_king = load_tf_seed_king_only()
    bridge = load_bridge()

    atlases = []
    for path, name in ((FINCHER_PATH, "fincher"), (PLASS_PATH, "plass"),
                       (CUI_PATH, "cui")):
        if not path.exists():
            print(f"  {name}: (missing, excluded from null)")
            continue
        adata = ad.read_h5ad(path)
        adata = prepare_atlas(adata, name, tf_ids_king)
        atlases.append((adata, name))
        print(f"  {name}: {adata.n_obs} cells x {adata.n_vars} genes, "
              f"leiden={adata.obs['leiden'].nunique()}")
    # Restrict each atlas to the TF gene panel (the pipeline's scoring
    # universe; keeps Wilcoxon cost bounded without changing the null).
    # --candidates mode: further restrict to the tested family's genes —
    # the Wilcoxon test on a gene is unaffected by which OTHER genes sit
    # in the panel (scanpy's rank_genes_groups is per-gene), so DE p's
    # for the tested genes are identical while the per-permutation cost
    # drops ~100x (the full-universe n=30 run took 8,614s; n=1000 at
    # family-restricted scale is tractable).
    if args.candidates:
        neural_csv = RUN_DIR / "rank_neural.csv"
        if not neural_csv.exists():
            print(f"--candidates requires {neural_csv}; run the pipeline first")
            return 1
        neural_ids = set(pd.read_csv(neural_csv)["gene_id"].astype(str))
        candidate_ids = [g for g in real_scores if g in neural_ids]
        print(f"--candidates mode: restricting p-value family to "
              f"{len(candidate_ids)} neural candidates "
              f"(atlas permutation unchanged)")
        tested_ids = set(candidate_ids)
        panel_target = tested_ids | {g + "_1" for g in tested_ids} | \
            {g.removesuffix("_1") for g in tested_ids}
    else:
        panel_target = tf_ids
        candidate_ids = list(real_scores.keys())
    panels = []
    for adata, name in atlases:
        if name == "fincher":
            v6_of = {v: bridge.v4_to_v6(v) for v in adata.var_names}
            keep = [v for v, v6 in v6_of.items() if v6 in panel_target]
        else:
            keep = [v for v in adata.var_names
                    if v in panel_target or (v + "_1") in panel_target]
        panel = adata[:, keep].copy()
        # CRITICAL performance fix (2026-09-06): prepare_atlas sets
        # adata.raw = adata (full atlas, 19-28k genes). scanpy's
        # rank_genes_groups and the log2FC helper both read .raw when
        # present, so every permutation silently ran Wilcoxon over the
        # FULL gene matrix despite the 100x smaller panel — the targeted
        # n=1000 run could not finish. The panel's DE inputs are
        # per-gene, so restricting raw to the panel changes nothing
        # statistically (identical p-values for the kept genes).
        panel.raw = panel
        panels.append((panel, name))
        print(f"  TF panel {name}: {len(keep)} genes")

    null_scores = {v6: [] for v6 in candidate_ids}
    for perm in range(args.n_perm):
        perm_expr = {}
        perm_spec = {}
        for apanel, name in panels:
            scores = permuted_atlas_streams(apanel, name, tf_ids, bridge, rng)
            for v6, (e, s) in scores.items():
                perm_expr[v6] = max(perm_expr.get(v6, 0.0), e)
                perm_spec[v6] = max(perm_spec.get(v6, 0.0), s)

        for v6 in candidate_ids:
            S = observed[v6].copy()
            ke = king_expr_floor.get(v6)
            ks = king_spec_floor.get(v6)
            pe = perm_expr.get(v6)
            ps = perm_spec.get(v6)

            # Expression: stream PRESENT iff permuted atlas hit OR King
            # floor (mirrors the pipeline); value = max(atlas, King).
            if pe is not None or ke is not None:
                S[0] = max(pe or 0.0, ke or 0.0)
            else:
                S[0] = np.nan
            # Specificity: same presence rule with the King breadth floor.
            if ps is not None or ks is not None:
                S[1] = max(ps or 0.0, ks or 0.0)
            else:
                S[1] = np.nan

            null_scores[v6].append(integrated_score_with_renorm(S, W_DEFAULT))

        if (perm + 1) % 5 == 0 or (perm + 1) == args.n_perm:
            print(f"  Completed {perm+1}/{args.n_perm} permutations")

    print("\n=== Empirical P-values (exchangeable-null, same statistic) ===")
    pvals = {}
    # family = candidate_ids (the 143-gene neural set in --candidates
    # mode; the full universe otherwise). out_df MUST be built from this
    # family, not real_scores — building from real_scores silently
    # exported all 11,696 genes with p=1.0 (no null draws) in --candidates
    # mode, poisoning the BH family.
    #
    # 2026-09-21 FIX: untestable genes (real <= ceiling + eps) get p=1.0
    # by enforcement — the >= counting for these genes was producing
    # floor p-values due to floating-point ties (e.g., surplus of 1e-16).
    for v6 in candidate_ids:
        real_s = real_scores[v6]
        is_untestable = real_s <= label_independent.get(v6, 0.0) + 1e-9
        if is_untestable:
            pvals[v6] = 1.0
        elif null_scores.get(v6):
            null_dist = np.array(null_scores[v6])
            pvals[v6] = (np.sum(null_dist >= real_s) + 1) / (len(null_dist) + 1)
        else:
            pvals[v6] = 1.0

    out_df = pd.DataFrame({
        "gene_id": list(pvals.keys()),
        "real_integrated_score": [real_scores.get(v, 0) for v in pvals],
        "empirical_p": list(pvals.values()),
        "n_perm": args.n_perm,
        "label_independent_ceiling": [label_independent.get(v, 0.0) for v in pvals],
        "untestable_by_permutation": [
            real_scores.get(v, 0) <= label_independent.get(v, 0.0) + 1e-9
            for v in pvals
        ],
    })
    real_names = dict(zip(real_rank["gene_id"], real_rank["gene_name"]))
    out_df["gene_name"] = out_df["gene_id"].map(real_names)
    out_df = out_df.sort_values(
        ["empirical_p", "real_integrated_score", "gene_id"],
        ascending=[True, False, True],
    )

    out_path = RESULTS_DIR / ("permutation_pvalues_neural.csv"
                              if args.candidates
                              else "permutation_pvalues_full.csv")

    # 2026-09-21 FIX: BH-FDR is applied ONLY over testable genes.
    # Untestable genes (p=1.0 by construction) would dilute the tested
    # family's q-values if included — 81 degenerate p=1.0 values plus
    # 12 float-tie floor p-values inflated the previous BH correction.
    # Untestable genes receive q=NaN (never tested) and
    # significant_fdr_05=False.
    from statsmodels.stats.multitest import multipletests
    testable_mask = ~out_df["untestable_by_permutation"].values
    out_df["q_bh_family"] = np.nan
    out_df["significant_fdr_05"] = False
    testable_pvals = out_df.loc[testable_mask, "empirical_p"].values
    if len(testable_pvals) > 0:
        _, qvals_testable, _, _ = multipletests(
            testable_pvals, alpha=0.05, method="fdr_bh"
        )
        out_df.loc[testable_mask, "q_bh_family"] = qvals_testable
        out_df.loc[testable_mask, "significant_fdr_05"] = qvals_testable < 0.05

    out_df.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")

    p_floor = 1.0 / (args.n_perm + 1)
    n_total = len(out_df)
    n_untestable = int(out_df["untestable_by_permutation"].sum())
    n_testable = n_total - n_untestable
    n_sig = int(out_df["significant_fdr_05"].sum())

    print(f"\nTestable genes: {n_testable} / {n_total} "
          f"({n_untestable} untestable — p=1.0 by construction, "
          f"q=NaN, excluded from BH family)")
    print(f"Significant after BH-FDR q<0.05 within the {n_testable}-gene "
          f"testable family: {n_sig} / {n_testable}")
    print(f"Resolution floor at n={args.n_perm}: min p = {p_floor:.5f} "
          f"(family BH q<0.05 {'reachable' if p_floor < 0.05 else 'NOT reachable'})")

    testable_df = out_df[~out_df["untestable_by_permutation"]]
    for _, row in testable_df.head(10).iterrows():
        nm = row["gene_name"]
        nm = nm if isinstance(nm, str) and nm and str(nm) != "nan" else row["gene_id"]
        print(f"  {str(nm)[:28]:>28} ({row['gene_id']}): "
              f"p={row['empirical_p']:.4f}, q={row['q_bh_family']:.4f}, "
              f"score={row['real_integrated_score']:.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
