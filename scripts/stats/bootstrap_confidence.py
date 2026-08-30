#!/usr/bin/env python
"""Bootstrap confidence intervals for integrated evidence scores.

Draws n=1000 bootstrap resamples from the evidence stream scores of all
candidates, computing the 95% CI on the integrated score for each gene.

Usage:
    python scripts/stats/bootstrap_confidence.py --n-boot 1000 --seed 42
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


def load_score_matrix():
    """Load the full candidate table and extract evidence stream columns."""
    candidates_path = RESULTS_DIR / "supplementary_table_S2_fixed_all249.csv"
    if not candidates_path.exists():
        candidates_path = RUN_DIR / "rank.csv"
    if not candidates_path.exists():
        raise FileNotFoundError("No candidate score file found in RESULTS_DIR or RUN_DIR")
    df = pd.read_csv(candidates_path)
    return df


def integrated_score(S, W):
    """Compute integrated score with missing-data renormalization."""
    mask = ~np.isnan(S)
    if not mask.any():
        return 0.0
    S_filled = np.where(np.isnan(S), 0.0, S)
    num = np.sum(S_filled * W)
    den = np.sum(W[mask])
    return num / den if den > 0 else 0.0


def bootstrap_ci(scores_matrix, weights, n_boot=1000, seed=42):
    """Return bootstrap means and 95% CIs for each candidate."""
    rng = np.random.default_rng(seed)
    n_candidates = scores_matrix.shape[0]
    boot_scores = np.zeros((n_boot, n_candidates))

    for b in range(n_boot):
        idx = rng.choice(n_candidates, size=n_candidates, replace=True)
        resampled = scores_matrix[idx]
        for j in range(n_candidates):
            boot_scores[b, j] = integrated_score(resampled[j], weights)

    means = np.mean(boot_scores, axis=0)
    ci_lo = np.percentile(boot_scores, 2.5, axis=0)
    ci_hi = np.percentile(boot_scores, 97.5, axis=0)
    return means, ci_lo, ci_hi


def main():
    parser = argparse.ArgumentParser(description="Bootstrap 95% CI for integrated scores")
    parser.add_argument("--n-boot", type=int, default=1000, help="Number of bootstrap resamples")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    print(f"=== Bootstrap Confidence Intervals (n={args.n_boot}) ===")

    df = load_score_matrix()
    stream_cols = [c for c in STREAMS if c in df.columns]
    if len(stream_cols) < 4:
        print(f"Warning: only {len(stream_cols)} stream columns found: {stream_cols}")

    scores_matrix = df[stream_cols].values.astype(float)
    weights = W_DEFAULT[:len(stream_cols)]
    weights = weights / weights.sum()

    print(f"Candidates: {scores_matrix.shape[0]}, Streams: {scores_matrix.shape[1]}")
    print(f"Running {args.n_boot} bootstrap resamples...")

    means, ci_lo, ci_hi = bootstrap_ci(scores_matrix, weights, args.n_boot, args.seed)

    out = df[["gene_id", "gene_name"]].copy() if "gene_name" in df.columns else df[["gene_id"]].copy()
    out["bootstrap_mean"] = means
    out["mean_score"] = means
    out["ci_95_lo"] = ci_lo
    out["ci_low"] = ci_lo
    out["ci_95_hi"] = ci_hi
    out["ci_high"] = ci_hi
    out["ci_width"] = ci_hi - ci_lo
    out = out.sort_values("bootstrap_mean", ascending=False)

    out_path = RESULTS_DIR / "bootstrap_scores_ci.csv"
    out.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")

    top10 = out.head(10)
    fig, ax = plt.subplots(figsize=(10, 6))
    y_pos = np.arange(len(top10))
    ax.barh(y_pos, top10["bootstrap_mean"], xerr=[
        top10["bootstrap_mean"] - top10["ci_95_lo"],
        top10["ci_95_hi"] - top10["bootstrap_mean"]
    ], capsize=3, color="#4C72B0", edgecolor="black", linewidth=0.5)
    ax.set_yticks(y_pos)
    label_col = "gene_name" if "gene_name" in top10.columns else "gene_id"
    ax.set_yticklabels(top10[label_col].values)
    ax.set_xlabel("Bootstrap Mean Integrated Score")
    ax.set_title(f"Top-10 Candidates: 95% Bootstrap CI (n={args.n_boot})")
    ax.invert_yaxis()
    plt.tight_layout()
    fig_path = FIG_DIR / "bootstrap_ci_top10.png"
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)
    print(f"Saved: {fig_path}")

    print(f"\nTop-10 by bootstrap mean:")
    for _, row in top10.iterrows():
        print(f"  {row.get('gene_name', row['gene_id']):>12}  "
              f"mean={row['bootstrap_mean']:.4f}  "
              f"95% CI=[{row['ci_95_lo']:.4f}, {row['ci_95_hi']:.4f}]")

    return 0


if __name__ == "__main__":
    sys.exit(main())
