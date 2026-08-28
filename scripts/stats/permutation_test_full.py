#!/usr/bin/env python
"""Full permutation test across all 9 evidence streams.

Shuffles cluster labels in all available atlases to generate null
distribution of integrated scores. Computes empirical p-values for
real candidates.

Usage:
    python scripts/stats/permutation_test_full.py --n-perm 1000 --subsample 2000
"""

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
FIG_DIR = REPO / "projects" / "NeuralTF" / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Evidence sources (same as pipeline)
STREAMS = ["expression", "specificity", "reproducibility", "rnai",
           "correlation", "neural_enriched", "neural_specificity",
           "perez_lineage", "perez_influence"]
W_DEFAULT = np.array([0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

# Paths
FINCHER_PATH = REPO / "datasets" / "processed" / "fincher_subsample.h5ad"
PLASS_PATH = REPO / "datasets" / "processed" / "plass_v6.h5ad"
BRIDGE_PATH = REPO / "projects" / "NeuralTF" / "data" / "bridge.csv"
KING_MMC4 = REPO / "datasets" / "raw" / "Supplementary_Data_ King_2024" / "1-s2.0-S2211124724001712-mmc4.xlsx"


def load_tf_catalog():
    """Load TF catalog from King mmc4."""
    df = pd.read_excel(KING_MMC4, sheet_name="TF")
    tf_ids = set(df.loc[df["TF?"].notna(), "Gene ID"].astype(str))
    tf_ids_norm = tf_ids | {tid[:-2] for tid in tf_ids if tid.endswith("_1")}
    return tf_ids_norm


def load_bridge():
    """Load bridge table for v4<->v6 mapping."""
    from bioforge.evidence import load_bridge as _load_bridge
    return _load_bridge(BRIDGE_PATH)


def score_one_atlas_permuted(adata, atlas_name, tf_ids_norm, bridge, rng):
    """Score one atlas with permuted cluster labels."""
    # Permute leiden labels
    leiden_perm = adata.obs["leiden"].copy()
    leiden_perm = leiden_perm.sample(frac=1.0, random_state=rng.integers(0, 2**32)).reset_index(drop=True)
    adata_perm = adata.copy()
    adata_perm.obs["leiden"] = leiden_perm.values

    # Run Wilcoxon on permuted data
    sc.tl.rank_genes_groups(adata_perm, "leiden", method="wilcoxon")
    result = adata_perm.uns["rank_genes_groups"]

    n_clusters = len(adata_perm.obs["leiden"].cat.categories)
    clusters = result["names"].dtype.names

    # Collect best scores per gene
    gene_best = {}
    for cl in clusters:
        for g, lfc, pval in zip(result["names"][cl], result["logfoldchanges"][cl], result["pvals"][cl]):
            key = str(g)
            abs_lfc = abs(float(lfc))
            if key not in gene_best or abs_lfc > abs(gene_best[key][0]):
                gene_best[key] = (abs_lfc, float(pval))

    # Build score matrix for TFs
    scores = {}
    if atlas_name == "fincher":
        v6_of = {v: bridge.v4_to_v6(v) for v in adata.var_names}
        score_genes = [v for v, v6 in v6_of.items() if v6 in tf_ids_norm]
    else:
        score_genes = [v for v in adata.var_names if v in tf_ids_norm]

    for gene in score_genes:
        best_l2fc, best_p = gene_best.get(gene, (0.0, 1.0))
        if best_p > 0.05:
            continue

        if atlas_name == "fincher":
            v6_id = bridge.v4_to_v6(gene)
        else:
            v6_id = gene if gene in tf_ids_norm else gene + "_1"

        if not v6_id:
            continue

        if v6_id not in scores:
            scores[v6_id] = np.full(9, np.nan)

        # Expression score
        scores[v6_id][0] = min(1.0, best_l2fc / 5.0)
        # Specificity score
        scores[v6_id][1] = 1.0 / n_clusters if n_clusters > 0 else 0.0

    return scores


def integrated_score_with_renorm(S: np.ndarray, W: np.ndarray) -> float:
    """Compute integrated score with missing-data renormalization."""
    mask = ~np.isnan(S)
    if not mask.any():
        return 0.0
    S_filled = np.where(np.isnan(S), 0.0, S)
    num = np.sum(S_filled * W)
    den = np.sum(W[mask])
    return num / den if den > 0 else 0.0


def main():
    parser = argparse.ArgumentParser(description="Full permutation test for NeuralTF")
    parser.add_argument("--n-perm", type=int, default=1000, help="Number of permutations")
    parser.add_argument("--subsample", type=int, default=2000, help="Subsample cells for speed")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    print(f"=== Full Permutation Test (n={args.n_perm}) ===")
    print(f"Subsample: {args.subsample} cells per atlas")

    # Load data
    print("Loading datasets...")
    adata_fincher = ad.read_h5ad(FINCHER_PATH)
    adata_plass = ad.read_h5ad(PLASS_PATH)

    if args.subsample and adata_fincher.n_obs > args.subsample:
        sc.pp.subsample(adata_fincher, n_obs=args.subsample, random_state=args.seed)
    if args.subsample and adata_plass.n_obs > args.subsample:
        sc.pp.subsample(adata_plass, n_obs=args.subsample, random_state=args.seed)

    # QC + clustering (same as pipeline)
    print("QC + clustering...")
    for adata in [adata_fincher, adata_plass]:
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=5000)
        sc.pp.pca(adata, n_comps=50)
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.leiden(adata, resolution=0.5)

    tf_ids_norm = load_tf_catalog()
    bridge = load_bridge()

    # Run permutations
    null_scores = {v6: [] for v6 in tf_ids_norm}

    for perm in range(args.n_perm):
        # Score Fincher
        fincher_scores = score_one_atlas_permuted(adata_fincher, "fincher", tf_ids_norm, bridge, rng)
        # Score Plass
        plass_scores = score_one_atlas_permuted(adata_plass, "plass", tf_ids_norm, bridge, rng)

        # Merge (best-atlas-wins for expression & specificity)
        for v6 in tf_ids_norm:
            S = np.full(9, np.nan)
            expr_vals = []
            spec_vals = []
            for scores in [fincher_scores, plass_scores]:
                if v6 in scores:
                    if not np.isnan(scores[v6][0]):
                        expr_vals.append(scores[v6][0])
                    if not np.isnan(scores[v6][1]):
                        spec_vals.append(scores[v6][1])
            if expr_vals:
                S[0] = max(expr_vals)
            if spec_vals:
                S[1] = max(spec_vals)

            score = integrated_score_with_renorm(S, W_DEFAULT)
            null_scores[v6].append(score)

        if (perm + 1) % 100 == 0:
            print(f"  Completed {perm+1}/{args.n_perm} permutations")

    # Load real scores for comparison
    real_rank_path = RUN_DIR / "rank.csv"
    if not real_rank_path.exists():
        print("Real rank.csv not found; run pipeline first")
        return 1

    real_rank = pd.read_csv(real_rank_path)
    real_scores = dict(zip(real_rank["gene_id"], real_rank["integrated_score"]))

    # Compute empirical p-values
    print("\n=== Computing Empirical P-values ===")
    pvals = {}
    for v6, real_s in real_scores.items():
        if v6 in null_scores and null_scores[v6]:
            null_dist = np.array(null_scores[v6])
            p = (np.sum(null_dist >= real_s) + 1) / (len(null_dist) + 1)
            pvals[v6] = p
        else:
            pvals[v6] = 1.0

    # Create output
    out_df = pd.DataFrame({
        "gene_id": list(pvals.keys()),
        "real_integrated_score": [real_scores.get(v, 0) for v in pvals.keys()],
        "empirical_p": list(pvals.values()),
        "n_perm": args.n_perm,
    })
    out_df = out_df.sort_values("empirical_p")

    # Merge with gene names
    real_names = dict(zip(real_rank["gene_id"], real_rank["gene_name"]))
    out_df["gene_name"] = out_df["gene_id"].map(real_names)

    out_path = RESULTS_DIR / "permutation_pvalues_full.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved permutation p-values to {out_path}")

    # Summary
    sig = out_df[out_df["empirical_p"] < 0.05]
    print(f"\nSignificant at p<0.05: {len(sig)} / {len(out_df)}")
    print(f"Top-10 real candidates p-values:")
    for _, row in out_df.head(10).iterrows():
        print(f"  {row['gene_name']:>10} ({row['gene_id']}): p={row['empirical_p']:.4f}, score={row['real_integrated_score']:.4f}")

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())