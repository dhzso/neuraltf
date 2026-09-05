#!/usr/bin/env python
"""Permutation test on the SAME multi-stream statistic the pipeline ranks.

2026-09-04 FULL REDESIGN:

The previous null was degenerate: it permuted only Fincher/Plass labels,
silently RETAINED the observed (King/Cui-inflated) expression for genes
with no permuted DE hit, forced null specificity to 1.0 via a dead regex,
and ran on 2,000-cell TF-only re-clustered atlases. 99.6% of genes landed
at the floor p-value - the test measured nothing.

Exchangeable-null design (this version):
- Permute cluster labels in ALL THREE scRNA atlases (Fincher, Plass, Cui)
  at PRODUCTION scale, re-running the pipeline's exact QC -> HVG+forced
  King TFs -> PCA -> neighbors -> igraph-leiden clustering, then Wilcoxon
  DE with the same global BH (q<=0.1) gate over genes x clusters.
- Label-DEPENDENT evidence (atlas DE expression + atlas specificity) is
  recomputed from the permuted structure alone.
- Label-INDEPENDENT evidence (King mmc7 table floors, RNAi mmc5, TF-pair
  correlations mmc6, Perez MOESM5/19 classes, neural gates from the King
  table) is held FIXED in both the real and null statistics - permuting
  cluster labels cannot change an external table.
- The King expression floor min(1, max_log2FC/5) and King fractional
  breadth specificity floor 1-(n-1)/(N-1) enter BOTH real and null via
  max(atlas value, King floor), exactly as the pipeline integrates them.
  Consequently a gene whose score is fully explained by table evidence
  (e.g. King log2FC >= 5 saturating expression at 1.0) has null == real
  by construction: its p-value is 1.0 and it is flagged
  label_independent_score=True - the honest statement that its rank
  rests on external tables, not on cluster-specific atlas expression.
- Stream PRESENCE mirrors the pipeline: expression/specificity are
  present in the null iff the gene has permuted atlas evidence OR a King
  floor; otherwise the streams are absent (NaN) as in the real score.
- Empirical p = (#{null >= real} + 1) / (n + 1), one-tailed add-one.

The null answers: "could this gene's integrated score arise from random
cluster structure in the atlases, given its label-independent evidence?"

Usage:
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
           "perez_lineage", "perez_influence"]
W_DEFAULT = np.array([0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

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
    """Replicate the pipeline's run_qc exactly (filter, normalize, log1p,
    HVG+forced King TFs, PCA, neighbors, igraph leiden 0.5, random_state
    fixed) so the permuted labels describe the SAME cluster structure the
    real scores were computed on."""
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
        S = np.full(9, np.nan)
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
    panels = []
    for adata, name in atlases:
        if name == "fincher":
            v6_of = {v: bridge.v4_to_v6(v) for v in adata.var_names}
            keep = [v for v, v6 in v6_of.items() if v6 in tf_ids]
        else:
            keep = [v for v in adata.var_names
                    if v in tf_ids or (v + "_1") in tf_ids]
        panels.append((adata[:, keep].copy(), name))
        print(f"  TF panel {name}: {len(keep)} genes")

    candidate_ids = list(real_scores.keys())

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
    for v6, real_s in real_scores.items():
        if null_scores.get(v6):
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

    out_path = RESULTS_DIR / "permutation_pvalues_full.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")

    testable = out_df[~out_df["untestable_by_permutation"]]
    sig = testable[testable["empirical_p"] < 0.05]
    print(f"\nSignificant at p<0.05 among testable genes: "
          f"{len(sig)} / {len(testable)} "
          f"({int(out_df['untestable_by_permutation'].sum())} untestable, "
          f"p=1.0 by construction)")
    for _, row in testable.head(10).iterrows():
        nm = row["gene_name"]
        nm = nm if isinstance(nm, str) and nm and str(nm) != "nan" else row["gene_id"]
        print(f"  {str(nm)[:28]:>28} ({row['gene_id']}): "
              f"p={row['empirical_p']:.4f}, score={row['real_integrated_score']:.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
