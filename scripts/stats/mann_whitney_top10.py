#!/usr/bin/env python
"""Mann-Whitney U test for top-10 vs remaining candidates.

Compares integrated scores of the top-10 ranked candidates against the
remaining candidates using the Mann-Whitney U test and rank-biserial
correlation effect size.

Usage:
    python scripts/stats/mann_whitney_top10.py
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


def rank_biserial_correlation(u_stat, n1, n2):
    """Compute rank-biserial correlation from Mann-Whitney U statistic."""
    return 1.0 - (2.0 * u_stat) / (n1 * n2)


def main():
    print("=== Mann-Whitney U Test: Top-10 vs Rest ===")

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

    n1, n2 = len(top10_scores), len(rest_scores)
    mean1, mean2 = np.mean(top10_scores), np.mean(rest_scores)
    median1, median2 = np.median(top10_scores), np.median(rest_scores)

    print(f"Top-10: n={n1}, mean={mean1:.4f}, median={median1:.4f}")
    print(f"Rest:   n={n2}, mean={mean2:.4f}, median={median2:.4f}")

    u_stat, p_value = stats.mannwhitneyu(top10_scores, rest_scores, alternative="greater")
    r_rb = rank_biserial_correlation(u_stat, n1, n2)

    u_two_sided, p_two_sided = stats.mannwhitneyu(top10_scores, rest_scores, alternative="two-sided")

    rank_sum_top10 = np.sum(stats.rankdata(np.concatenate([top10_scores, rest_scores]))[:n1])
    expected_rank_sum = n1 * (n1 + n2 + 1) / 2

    cliff_d = 0.0
    for t in top10_scores:
        cliff_d += np.sum(t > rest_scores) - np.sum(t < rest_scores)
    cliff_d /= (n1 * n2)

    results = {
        "test": "Mann-Whitney U (one-sided, greater)",
        "top10_n": int(n1),
        "rest_n": int(n2),
        "top10_mean": float(mean1),
        "rest_mean": float(mean2),
        "top10_median": float(median1),
        "rest_median": float(median2),
        "u_statistic": float(u_stat),
        "p_value_one_sided": float(p_value),
        "u_statistic_two_sided": float(u_two_sided),
        "p_value_two_sided": float(p_two_sided),
        "rank_biserial_correlation": float(r_rb),
        "cliffs_delta": float(cliff_d),
        "rank_sum_top10": float(rank_sum_top10),
        "expected_rank_sum_under_null": float(expected_rank_sum),
        "top10_genes": sorted(list(top10_genes)),
    }

    print(f"\n--- Results ---")
    print(f"  U statistic:           {u_stat:.1f}")
    print(f"  p-value (one-sided):   {p_value:.4e}")
    print(f"  p-value (two-sided):   {p_two_sided:.4e}")
    print(f"  Rank-biserial r:       {r_rb:.4f}")
    print(f"  Cliff's delta:         {cliff_d:.4f}")
    print(f"  Rank sum (top-10):     {rank_sum_top10:.0f}")
    print(f"  Expected under H0:     {expected_rank_sum:.0f}")

    print(f"\nTop-10 genes:")
    top10_df = df_sorted.head(10)
    for _, row in top10_df.iterrows():
        print(f"  {row.get('gene_name', row[gene_col]):>12s}  score={row[score_col]:.4f}")

    out_path = RESULTS_DIR / "mann_whitney_top10.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
