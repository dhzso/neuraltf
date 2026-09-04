#!/usr/bin/env python
"""Leave-one-atlas-out sensitivity analysis.

Re-ranks candidates with each evidence atlas removed, tracking top-10
stability across leave-one-out iterations.

Usage:
    python scripts/stats/leave_one_atlas_out.py
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
W_DEFAULT = np.array([0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])


def integrated_score(S, W):
    """Compute integrated score with missing-data renormalization."""
    mask = ~np.isnan(S)
    if not mask.any():
        return 0.0
    S_filled = np.where(np.isnan(S), 0.0, S)
    w_masked = W[mask]
    return np.sum(S_filled[mask] * w_masked) / w_masked.sum()


def main():
    print("=== Leave-One-Atlas-Out Analysis ===")

    p = RUN_DIR / "rank.csv"
    if not p.exists():
        print("Error: run the pipeline first (rank.csv missing)")
        return 1

    df = pd.read_csv(p).drop_duplicates(subset="gene_id", keep="first")

    gene_col = "gene_id" if "gene_id" in df.columns else "gene_id_v6"
    name_col = "gene_name" if "gene_name" in df.columns else gene_col
    score_col = "integrated_score" if "integrated_score" in df.columns else "composite_score"

    stream_cols = [s for s in STREAMS if s in df.columns]
    if len(stream_cols) < 3:
        stream_cols = [c for c in df.columns if c in STREAMS]

    scores_matrix = df[stream_cols].values.astype(float)
    weights = W_DEFAULT[:len(stream_cols)]
    weights = weights / weights.sum()

    baseline_scores = np.array([integrated_score(row, weights) for row in scores_matrix])
    baseline_order = np.argsort(-baseline_scores)
    top10_indices = baseline_order[:10]
    top10_ids = df[gene_col].values[top10_indices]
    top10_names = df[name_col].fillna(df[gene_col]).values[top10_indices]

    print(f"Candidates: {len(df)}, Streams: {len(stream_cols)}")
    print(f"Baseline top-10: {list(top10_names)}")

    results = []
    # Rank matrix for heatmap: rows = top10 candidates, cols = excluded streams
    rank_matrix_dict = {
        "gene_id": list(top10_ids),
        "gene_name": list(top10_names),
        "full_rank": list(range(1, 11)),
        "full_score": [float(baseline_scores[idx]) for idx in top10_indices],
    }

    for i, stream in enumerate(stream_cols):
        loo_streams = [s for j, s in enumerate(stream_cols) if j != i]
        loo_idx = [j for j in range(len(stream_cols)) if j != i]
        loo_weights = np.delete(weights, i)
        loo_weights = loo_weights / loo_weights.sum()
        loo_matrix = scores_matrix[:, loo_idx]

        loo_scores = np.array([integrated_score(row, loo_weights) for row in loo_matrix])
        loo_ranks = pd.Series(loo_scores).rank(ascending=False, method="min").values
        loo_order = np.argsort(-loo_scores)
        loo_top10 = set(df[gene_col].values[loo_order[:10]])

        overlap = set(top10_ids) & loo_top10
        jaccard = len(overlap) / len(set(top10_ids) | loo_top10) if len(set(top10_ids) | loo_top10) > 0 else 0
        rank_corr = pd.Series(baseline_scores).corr(pd.Series(loo_scores), method="spearman")

        # Record ranks of baseline top-10 when this stream is excluded
        rank_matrix_dict[stream] = [int(loo_ranks[idx]) for idx in top10_indices]

        results.append({
            "excluded_atlas": stream,
            "n_remaining_streams": len(loo_streams),
            "top10_overlap": len(overlap),
            "top10_jaccard": float(jaccard),
            "spearman_correlation": float(rank_corr),
            "overlap_genes": sorted(list(overlap)),
            "new_in_top10": sorted(list(loo_top10 - set(top10_ids))),
            "dropped_from_top10": sorted(list(set(top10_ids) - loo_top10)),
        })
        print(f"  Exclude {stream:>20s}: overlap={len(overlap)}/10, "
              f"Jaccard={jaccard:.3f}, rho={rank_corr:.4f}")

    # Write rank matrix to loo_atlas_stability.csv for figure 27
    matrix_df = pd.DataFrame(rank_matrix_dict)
    out_path = RESULTS_DIR / "loo_atlas_stability.csv"
    matrix_df.to_csv(out_path, index=False)
    print(f"\nSaved stability matrix: {out_path}")

    # Write summary metrics to loo_atlas_summary.csv
    summary_df = pd.DataFrame(results)
    summary_path = RESULTS_DIR / "loo_atlas_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"Saved summary metrics: {summary_path}")

    avg_overlap = np.mean([r["top10_overlap"] for r in results])
    avg_jaccard = np.mean([r["top10_jaccard"] for r in results])
    avg_rho = np.mean([r["spearman_correlation"] for r in results])
    print(f"\nMean overlap: {avg_overlap:.1f}/10")
    print(f"Mean Jaccard: {avg_jaccard:.3f}")
    print(f"Mean Spearman rho: {avg_rho:.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
