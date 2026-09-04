#!/usr/bin/env python
"""Permutation test on the SAME multi-stream statistic the pipeline ranks.

WS3 fix: the previous null model filled only expression/specificity from
permuted cluster labels (denominator 0.3 of the weight mass) while real
scores in rank.csv use all 9 streams including King/Perez/RNAi evidence
(denominator 1.0). Comparing real vs null across different score spaces
made 263/278 "significant" for reasons unrelated to DE signal.

Correct approach: keep the NON-DE-stream evidence (rnai, correlation,
neural_enriched, neural_specificity, perez_lineage, perez_influence,
reproducibility) FIXED per gene as observed, and permute only the
atlas-derived expression/specificity signal. The null then answers the
right question: "could this gene's integrated score arise from random
cluster structure, given its other evidence?"

The empirical p-value uses (count+1)/(n+1) — the add-one correction that
keeps p>0 and is unbiased under the null.

Usage:
    python scripts/stats/permutation_test_full.py --n-perm 100 --subsample 2000
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
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

STREAMS = ["expression", "specificity", "reproducibility", "rnai",
           "correlation", "neural_enriched", "neural_specificity",
           "perez_lineage", "perez_influence"]
W_DEFAULT = np.array([0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

FINCHER_PATH = REPO / "datasets" / "processed" / "fincher_subsample.h5ad"
PLASS_PATH = REPO / "datasets" / "processed" / "plass_v6.h5ad"
BRIDGE_PATH = REPO / "projects" / "NeuralTF" / "data" / "bridge.csv"
KING_MMC4 = (REPO / "datasets" / "raw" / "Supplementary_Data_ King_2024"
             / "1-s2.0-S2211124724001712-mmc4.xlsx")


def load_tf_catalog() -> set[str]:
    df = pd.read_excel(KING_MMC4, sheet_name="TF")
    tf_ids = set(df.loc[df["TF?"].notna(), "Gene ID"].astype(str))
    return tf_ids | {tid[:-2] for tid in tf_ids if tid.endswith("_1")}


def load_bridge():
    from bioforge.evidence import load_bridge as _load_bridge
    return _load_bridge(BRIDGE_PATH)


def score_one_atlas_permuted(adata, atlas_name, tf_ids_norm, bridge, rng,
                             n_sig_ref: dict[str, int]):
    """Score one atlas with permuted cluster labels on pre-filtered TF genes.

    Mirrors the pipeline: one-tailed max positive lfc among significant
    clusters (BH q<=0.1 via per-cluster rank p) for expression; specificity
    uses the observed significant-cluster count (permutation cannot create
    real breadth, so the observed n_sig is the honest denominator).
    """
    perm_labels = rng.permutation(adata.obs["leiden"].astype(str).values)
    adata.obs["leiden_perm"] = pd.Categorical(perm_labels)

    sc.tl.rank_genes_groups(adata, "leiden_perm", method="wilcoxon")
    result = adata.uns["rank_genes_groups"]
    clusters = result["names"].dtype.names

    gene_best = {}
    for cl in clusters:
        for g, lfc, pval in zip(
            result["names"][cl], result["logfoldchanges"][cl], result["pvals"][cl]
        ):
            key = str(g)
            pos_lfc = max(0.0, float(lfc))
            p = float(pval)
            # significance gate comparable to the pipeline's BH q<=0.1
            # (rank_genes_groups already reports raw p; the pipeline's
            # global BH over genes x clusters is approximated here by the
            # per-gene min-p since we filter to the TF subset)
            if p <= 0.05 and pos_lfc > 0:
                if key not in gene_best or pos_lfc > gene_best[key][0]:
                    gene_best[key] = (pos_lfc, p)

    scores = {}
    for gene in adata.var_names:
        best_l2fc, _ = gene_best.get(gene, (0.0, 1.0))
        if best_l2fc <= 0:
            continue
        if atlas_name == "fincher":
            v6_id = bridge.v4_to_v6(gene)
        else:
            v6_id = gene if gene in tf_ids_norm else gene + "_1"

        if not v6_id:
            continue
        # expression: one-tailed pseudo-lfc scaled as in the null's
        # comparison space (the DE score the pipeline would assign);
        # specificity: observed significant breadth, same as pipeline
        n_sig = max(n_sig_ref.get(v6_id, 1), 1)
        scores[v6_id] = (min(1.0, best_l2fc / 5.0), 1.0 / n_sig)
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
    parser = argparse.ArgumentParser(
        description="Permutation test on the full multi-stream statistic"
    )
    parser.add_argument("--n-perm", type=int, default=100)
    parser.add_argument("--subsample", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    print(f"=== Permutation Test on full statistic (n={args.n_perm}) ===")

    # Real scores + observed stream matrix (fixed under permutation)
    real_rank_path = RUN_DIR / "rank.csv"
    if not real_rank_path.exists():
        print("Real rank.csv not found; run pipeline first")
        return 1
    real_rank = pd.read_csv(real_rank_path).drop_duplicates(
        subset="gene_id", keep="first")
    real_scores = dict(zip(real_rank["gene_id"], real_rank["integrated_score"]))

    # Observed per-gene stream vectors (non-DE streams stay FIXED)
    stream_cols = [c for c in STREAMS if c in real_rank.columns]
    stream_idx = {c: STREAMS.index(c) for c in stream_cols}
    observed = {}
    n_sig_ref = {}
    for _, row in real_rank.iterrows():
        S = np.full(9, np.nan)
        for c in stream_cols:
            v = row[c]
            S[stream_idx[c]] = float(v) if pd.notna(v) else np.nan
        observed[row["gene_id"]] = S
        # observed significant-cluster count drives specificity denominator
        note = str(row.get("specificity", ""))
        try:
            import re as _re
            m = _re.search(r"n_sig_clusters=(\d+)", note)
            if m:
                n_sig_ref[row["gene_id"]] = int(m.group(1))
        except Exception:
            pass

    print(f"Real candidates: {len(real_rank)}")

    # Load atlases
    print("Loading atlases...")
    adata_fincher = ad.read_h5ad(FINCHER_PATH)
    adata_plass = ad.read_h5ad(PLASS_PATH)
    for adata in (adata_fincher, adata_plass):
        if args.subsample and adata.n_obs > args.subsample:
            sc.pp.subsample(adata, n_obs=args.subsample, random_state=args.seed)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=min(5000, adata.n_vars))
        sc.pp.pca(adata, n_comps=min(50, adata.n_vars - 1, adata.n_obs - 1))
        sc.pp.neighbors(adata, n_neighbors=10,
                        n_pcs=min(40, adata.n_vars - 1, adata.n_obs - 1))
        sc.tl.leiden(adata, resolution=0.5)

    tf_ids_norm = load_tf_catalog()
    bridge = load_bridge()

    # Pre-filter to TF genes for permutation speed
    v6_of_fincher = {v: bridge.v4_to_v6(v) for v in adata_fincher.var_names}
    fincher_tf_genes = [v for v, v6 in v6_of_fincher.items() if v6 in tf_ids_norm]
    adata_fincher_tf = adata_fincher[:, fincher_tf_genes].copy()
    plass_tf_genes = [v for v in adata_plass.var_names
                      if v in tf_ids_norm or (v + "_1") in tf_ids_norm]
    adata_plass_tf = adata_plass[:, plass_tf_genes].copy()
    print(f"TF-filtered: Fincher {adata_fincher_tf.n_vars}, Plass {adata_plass_tf.n_vars}")

    candidate_ids = list(real_scores.keys())

    # Permutations: permute ONLY expression/specificity; keep all other
    # observed streams fixed. Same 9-stream renormalized statistic.
    null_scores = {v6: [] for v6 in candidate_ids}
    for perm in range(args.n_perm):
        fincher_scores = score_one_atlas_permuted(
            adata_fincher_tf, "fincher", tf_ids_norm, bridge, rng, n_sig_ref)
        plass_scores = score_one_atlas_permuted(
            adata_plass_tf, "plass", tf_ids_norm, bridge, rng, n_sig_ref)

        for v6 in candidate_ids:
            S = observed.get(v6, np.full(9, np.nan)).copy()
            expr_vals, spec_vals = [], []
            for scores in (fincher_scores, plass_scores):
                if v6 in scores:
                    expr_vals.append(scores[v6][0])
                    spec_vals.append(scores[v6][1])
            if expr_vals:
                S[0] = max(expr_vals)   # best-atlas-wins, as the pipeline
            if spec_vals:
                S[1] = max(spec_vals)
            else:
                S[1] = np.nan if np.isnan(observed.get(v6, S)[1]) else 0.0
            null_scores[v6].append(integrated_score_with_renorm(S, W_DEFAULT))

        if (perm + 1) % 25 == 0 or (perm + 1) == args.n_perm:
            print(f"  Completed {perm+1}/{args.n_perm} permutations")

    # Empirical p-values with add-one correction
    print("\n=== Empirical P-values (same-statistic null) ===")
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
    })
    real_names = dict(zip(real_rank["gene_id"], real_rank["gene_name"]))
    out_df["gene_name"] = out_df["gene_id"].map(real_names)
    out_df = out_df.sort_values("empirical_p")

    out_path = RESULTS_DIR / "permutation_pvalues_full.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")

    sig = out_df[out_df["empirical_p"] < 0.05]
    print(f"\nSignificant at p<0.05: {len(sig)} / {len(out_df)}")
    for _, row in out_df.head(10).iterrows():
        print(f"  {str(row['gene_name']):>12} ({row['gene_id']}): "
              f"p={row['empirical_p']:.4f}, score={row['real_integrated_score']:.4f}")

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
