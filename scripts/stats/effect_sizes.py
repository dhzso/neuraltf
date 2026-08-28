#!/usr/bin/env python
"""Effect size analysis for score distributions.

Computes Cliff's delta for top-10 vs rest, Cohen's d for neural vs
non-neural, and Mann-Whitney U test for group separations.

Usage:
    python scripts/stats/effect_sizes.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
FIG_DIR = REPO / "projects" / "NeuralTF" / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


def cliffs_delta(x, y):
    """Compute Cliff's delta effect size. Returns delta in [-1, 1]."""
    n_x, n_y = len(x), len(y)
    if n_x == 0 or n_y == 0:
        return 0.0
    dominance = 0
    for xi in x:
        dominance += np.sum(xi > y) - np.sum(xi < y)
    return dominance / (n_x * n_y)


def cohens_d(x, y):
    """Compute Cohen's d pooled effect size."""
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return 0.0
    pooled_std = np.sqrt(((nx - 1) * np.var(x, ddof=1) + (ny - 1) * np.var(y, ddof=1)) / (nx + ny - 2))
    if pooled_std == 0:
        return 0.0
    return (np.mean(x) - np.mean(y)) / pooled_std


def main():
    print("=== Effect Size Analysis ===")

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

    score_col = None
    for c in ["integrated_score", "composite_score", "dirichlet_median_score", "fixed_weight_score"]:
        if c in df.columns:
            score_col = c
            break
    if score_col is None:
        print("Error: no score column found")
        return 1

    gene_col = "gene_id" if "gene_id" in df.columns else "gene_id_v6"
    df = df.dropna(subset=[score_col])

    df_sorted = df.sort_values(score_col, ascending=False)
    top10_genes = set(df_sorted[gene_col].head(10).values)
    top10_scores = df.loc[df[gene_col].isin(top10_genes), score_col].values
    rest_scores = df.loc[~df[gene_col].isin(top10_genes), score_col].values

    neural_mask = df["proof_status"] == "known_rnai_validated"
    neural_scores = df.loc[neural_mask, score_col].values
    non_neural_scores = df.loc[~neural_mask, score_col].values

    results = {}

    print(f"\n--- Top-10 vs Rest ---")
    print(f"  Top-10: n={len(top10_scores)}, mean={np.mean(top10_scores):.4f}, median={np.median(top10_scores):.4f}")
    print(f"  Rest:   n={len(rest_scores)}, mean={np.mean(rest_scores):.4f}, median={np.median(rest_scores):.4f}")

    cd_top10 = cliffs_delta(top10_scores, rest_scores)
    d_top10 = cohens_d(top10_scores, rest_scores)
    u_top10, p_top10 = stats.mannwhitneyu(top10_scores, rest_scores, alternative="greater")

    results["top10_vs_rest"] = {
        "cliffs_delta": float(cd_top10),
        "cohens_d": float(d_top10),
        "mann_whitney_u": float(u_top10),
        "p_value": float(p_top10),
        "top10_n": int(len(top10_scores)),
        "rest_n": int(len(rest_scores)),
        "top10_mean": float(np.mean(top10_scores)),
        "rest_mean": float(np.mean(rest_scores)),
    }
    print(f"  Cliff's delta: {cd_top10:.4f}")
    print(f"  Cohen's d:     {d_top10:.4f}")
    print(f"  Mann-Whitney U: {u_top10:.1f}, p={p_top10:.4e}")

    print(f"\n--- Neural vs Non-Neural ---")
    print(f"  Neural:     n={len(neural_scores)}, mean={np.mean(neural_scores):.4f}")
    print(f"  Non-neural: n={len(non_neural_scores)}, mean={np.mean(non_neural_scores):.4f}")

    cd_neural = cliffs_delta(neural_scores, non_neural_scores)
    d_neural = cohens_d(neural_scores, non_neural_scores)
    u_neural, p_neural = stats.mannwhitneyu(neural_scores, non_neural_scores, alternative="greater")

    results["neural_vs_non_neural"] = {
        "cliffs_delta": float(cd_neural),
        "cohens_d": float(d_neural),
        "mann_whitney_u": float(u_neural),
        "p_value": float(p_neural),
        "neural_n": int(len(neural_scores)),
        "non_neural_n": int(len(non_neural_scores)),
        "neural_mean": float(np.mean(neural_scores)),
        "non_neural_mean": float(np.mean(non_neural_scores)),
    }
    print(f"  Cliff's delta: {cd_neural:.4f}")
    print(f"  Cohen's d:     {d_neural:.4f}")
    print(f"  Mann-Whitney U: {u_neural:.1f}, p={p_neural:.4e}")

    out_path = RESULTS_DIR / "effect_sizes.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
