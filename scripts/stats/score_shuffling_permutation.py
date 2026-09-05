#!/usr/bin/env python
"""Score-shuffling permutation test (alternative null model).

Shuffles evidence stream assignments across candidates rather than
permuting cluster labels. This tests the null hypothesis that stream
assignments are exchangeable.

Usage:
    python scripts/stats/score_shuffling_permutation.py --n-perm 1000 --seed 42
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
FIG_DIR = REPO / "projects" / "NeuralTF" / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

STREAMS = [
    "expression", "specificity", "reproducibility", "rnai",
    "correlation", "neural_enriched", "neural_specificity",
    "perez_lineage", "perez_influence",
]
# Must match bioforge.evidence.scoring.DEFAULT_WEIGHTS exactly
# (the old perez_lineage=0.2 made "real" scores differ from rank.csv).
W_DEFAULT = np.array([0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])


def compute_all_integrated_scores(scores: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Vectorized calculation of integrated scores across all candidates."""
    scores_filled = np.nan_to_num(scores, nan=0.0)
    valid_mask = ~np.isnan(scores)
    w_denom = np.dot(valid_mask, weights)
    w_denom = np.where(w_denom > 0, w_denom, 1.0)
    return np.dot(scores_filled, weights) / w_denom


def shuffle_streams(scores_matrix, rng):
    """Shuffle each evidence stream independently across candidates."""
    shuffled = scores_matrix.copy()
    for j in range(shuffled.shape[1]):
        col = shuffled[:, j]
        valid = ~np.isnan(col)
        if valid.sum() > 1:
            shuffled[valid, j] = col[valid][rng.permutation(valid.sum())]
    return shuffled


def main():
    parser = argparse.ArgumentParser(description="Score-shuffling permutation test")
    parser.add_argument("--n-perm", type=int, default=200, help="Number of permutations (default: 200)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    print(f"=== Score-Shuffling Permutation Test (n={args.n_perm}) ===")

    rng = np.random.default_rng(args.seed)

    p = RUN_DIR / "rank.csv"
    if not p.exists():
        print("Error: no candidate score file found (run the pipeline first)")
        return 1

    df = pd.read_csv(p).drop_duplicates(subset="gene_id", keep="first")

    gene_col = "gene_id" if "gene_id" in df.columns else "gene_id_v6"
    score_col = "integrated_score" if "integrated_score" in df.columns else "composite_score"

    stream_cols = [s for s in STREAMS if s in df.columns]
    if len(stream_cols) < 3:
        stream_cols = [c for c in df.columns if c in STREAMS]

    scores_matrix = df[stream_cols].values.astype(float)
    weights = W_DEFAULT[:len(stream_cols)]
    weights = weights / weights.sum()

    real_scores = compute_all_integrated_scores(scores_matrix, weights)
    df["real_integrated_score"] = real_scores
    df_sorted = df.sort_values("real_integrated_score", ascending=False)

    print(f"Candidates: {len(df)}, Streams: {len(stream_cols)}")
    print(f"Top-10 real scores: {df_sorted['real_integrated_score'].head(10).values}")

    null_distributions = {gene: [] for gene in df[gene_col].values}
    top10_null_max = []

    for perm in range(args.n_perm):
        shuffled = shuffle_streams(scores_matrix, rng)
        perm_scores = compute_all_integrated_scores(shuffled, weights)
        for idx, gene in enumerate(df[gene_col].values):
            null_distributions[gene].append(perm_scores[idx])
        top10_null_max.append(np.max(perm_scores))

        if (perm + 1) % 50 == 0 or (perm + 1) == args.n_perm:
            print(f"  Completed {perm+1}/{args.n_perm} permutations")

    print("\n=== Computing Empirical P-values ===")
    pvals = []
    for _, row in df_sorted.iterrows():
        gene = row[gene_col]
        real_s = row["real_integrated_score"]
        null_dist = np.array(null_distributions[gene])
        p = (np.sum(null_dist >= real_s) + 1) / (args.n_perm + 1)
        pvals.append(p)

    df_out = df_sorted[[gene_col, "gene_name", "real_integrated_score"]].copy() if "gene_name" in df_sorted.columns else df_sorted[[gene_col, "real_integrated_score"]].copy()
    df_out = df_out.reset_index(drop=True)
    df_out["empirical_p_shuffled"] = pvals

    top10_real_max = np.max(real_scores)
    top10_p_global = (np.sum(np.array(top10_null_max) >= top10_real_max) + 1) / (args.n_perm + 1)

    out_path = RESULTS_DIR / "score_shuffling_pvalues.csv"
    df_out.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")

    n_sig_05 = sum(p < 0.05 for p in pvals)
    n_sig_01 = sum(p < 0.01 for p in pvals)
    print(f"\nSignificant at p<0.05: {n_sig_05}/{len(pvals)}")
    print(f"Significant at p<0.01: {n_sig_01}/{len(pvals)}")
    print(f"Global test (max real vs null max): p={top10_p_global:.4e}")

    print("\nTop-10 by score with shuffled-null p-values:")
    for i, (_, row) in enumerate(df_out.head(10).iterrows()):
        name = row.get("gene_name", row[gene_col])
        name = str(name) if (name is not None and str(name) != "nan" and str(name) != "None") else str(row[gene_col])
        print(f"  {i+1:>2d}. {name[:24]:>24}  score={row['real_integrated_score']:.4f}  "
              f"p={row['empirical_p_shuffled']:.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
