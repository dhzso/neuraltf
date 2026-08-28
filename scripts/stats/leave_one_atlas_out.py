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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
    "perez_lineage",
]
W_DEFAULT = np.array([0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2])


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

    for fname in ["supplementary_table_S2_fixed_all249.csv", "fstf_ranked_all.csv"]:
        p = RESULTS_DIR / fname
        if p.exists():
            df = pd.read_csv(p)
            break
    else:
        p = RUN_DIR / "rank.csv"
        if p.exists():
            df = pd.read_csv(p)
        else:
            print("Error: no candidate score file found")
            return 1

    gene_col = "gene_id" if "gene_id" in df.columns else "gene_id_v6"
    name_col = "gene_name" if "gene_name" in df.columns else None
    score_col = None
    for c in ["integrated_score", "composite_score"]:
        if c in df.columns:
            score_col = c
            break

    stream_cols = [s for s in STREAMS if s in df.columns]
    if len(stream_cols) < 3:
        print(f"Warning: only {len(stream_cols)} stream columns found: {stream_cols}")
        print("Using available numeric evidence columns instead")
        stream_cols = [c for c in df.columns if c in STREAMS or (
            df[c].dtype in [np.float64, np.float32, np.int64] and
            c not in [gene_col, name_col, score_col, "rank", "n_streams", "completeness"]
        )]

    scores_matrix = df[stream_cols].values.astype(float)
    weights = W_DEFAULT[:len(stream_cols)]
    weights = weights / weights.sum()

    baseline_scores = np.array([integrated_score(row, weights) for row in scores_matrix])
    baseline_order = np.argsort(-baseline_scores)
    baseline_top10 = set(df[gene_col].values[baseline_order[:10]])

    print(f"Candidates: {len(df)}, Streams: {len(stream_cols)}")
    print(f"Baseline top-10: {sorted(baseline_top10)}")

    results = []
    loo_top10_sets = {}

    for i, stream in enumerate(stream_cols):
        loo_streams = [s for j, s in enumerate(stream_cols) if j != i]
        loo_idx = [j for j in range(len(stream_cols)) if j != i]
        loo_weights = np.delete(weights, i)
        loo_weights = loo_weights / loo_weights.sum()
        loo_matrix = scores_matrix[:, loo_idx]

        loo_scores = np.array([integrated_score(row, loo_weights) for row in loo_matrix])
        loo_order = np.argsort(-loo_scores)
        loo_top10 = set(df[gene_col].values[loo_order[:10]])

        overlap = baseline_top10 & loo_top10
        jaccard = len(overlap) / len(baseline_top10 | loo_top10) if len(baseline_top10 | loo_top10) > 0 else 0
        rank_corr = pd.Series(baseline_scores).corr(pd.Series(loo_scores), method="spearman")

        loo_top10_sets[stream] = loo_top10

        results.append({
            "excluded_atlas": stream,
            "n_remaining_streams": len(loo_streams),
            "top10_overlap": len(overlap),
            "top10_jaccard": float(jaccard),
            "spearman_correlation": float(rank_corr),
            "overlap_genes": sorted(list(overlap)),
            "new_in_top10": sorted(list(loo_top10 - baseline_top10)),
            "dropped_from_top10": sorted(list(baseline_top10 - loo_top10)),
        })
        print(f"  Exclude {stream:>20s}: overlap={len(overlap)}/10, "
              f"Jaccard={jaccard:.3f}, rho={rank_corr:.4f}")

    out_df = pd.DataFrame(results)
    out_path = RESULTS_DIR / "loo_atlas_stability.csv"
    out_df.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")

    avg_overlap = np.mean([r["top10_overlap"] for r in results])
    avg_jaccard = np.mean([r["top10_jaccard"] for r in results])
    avg_rho = np.mean([r["spearman_correlation"] for r in results])
    print(f"\nMean overlap: {avg_overlap:.1f}/10")
    print(f"Mean Jaccard: {avg_jaccard:.3f}")
    print(f"Mean Spearman rho: {avg_rho:.4f}")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    excluded = [r["excluded_atlas"] for r in results]

    ax = axes[0]
    overlaps = [r["top10_overlap"] for r in results]
    ax.barh(excluded, overlaps, color="#4C72B0", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Top-10 Overlap with Baseline")
    ax.set_title("Stability: Top-10 Overlap")
    ax.set_xlim([0, 10])
    ax.axvline(x=avg_overlap, color="red", linestyle="--", label=f"Mean={avg_overlap:.1f}")
    ax.legend()

    ax = axes[1]
    jaccards = [r["top10_jaccard"] for r in results]
    ax.barh(excluded, jaccards, color="#55A868", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Jaccard Index")
    ax.set_title("Stability: Jaccard Similarity")
    ax.set_xlim([0, 1])

    ax = axes[2]
    rhos = [r["spearman_correlation"] for r in results]
    ax.barh(excluded, rhos, color="#C44E52", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Spearman rho")
    ax.set_title("Stability: Rank Correlation")
    ax.set_xlim([0.5, 1.0])

    plt.tight_layout()
    fig_path = FIG_DIR / "loo_atlas_stability.png"
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)
    print(f"Saved: {fig_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
