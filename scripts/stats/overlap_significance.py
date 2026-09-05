#!/usr/bin/env python
"""Overlap significance tests between ranking methods.

Tests whether the top-10 overlap between fixed, centered, and uniform
Dirichlet methods is greater than expected by chance using:
- Hypergeometric test for pairwise overlaps (primary statistic)
- Binomial pooled-rate test at the correct null p0 = 10/N
- Three-way overlap with a proper three-set null (NOT the pairwise null
  the previous version used, which was conservative by ~8 orders of
  magnitude for small overlaps)

2026-09-04 alignment: the "fixed" arm now loads the published fixed-method
shortlist (top10_neural_tfs_prioritized.csv) so all three arms are the
same quantity class (composite dual-track shortlists), matching
cross_method_correction.py.

Usage:
    python scripts/stats/overlap_significance.py
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


def hypergeometric_test(k, n, K, N):
    """P-value for observing >= k successes in sample of n from population of N with K successes."""
    p = stats.hypergeom.sf(k - 1, N, K, n)
    return p


def three_way_null_p(k3, n, N, n_methods=3):
    """Exact null p for k3 genes shared by all n_methods independent
    random n-subsets of an N-gene universe.

    P(a given gene is in one random n-subset) = n/N; independence across
    methods gives P(in all n_methods) = (n/N)^n_methods, and the count
    of 3-way-shared genes is Binomial(N, (n/N)^n_methods) under the
    independent-random-subsets null. One-sided P(X >= k3).
    """
    p_each = (n / N) ** n_methods
    if hasattr(stats, "binomtest"):
        return float(stats.binomtest(k3, N, p_each, alternative="greater").pvalue)
    return float(stats.binom_test(k3, N, p_each, alternative="greater"))


def _binom_p(k, n, p):
    """Compute binomial test p-value compatible with SciPy <1.12 and >=1.12."""
    if hasattr(stats, "binomtest"):
        return float(stats.binomtest(k, n, p, alternative="greater").pvalue)
    return float(stats.binom_test(k, n, p, alternative="greater"))


def main():
    print("=== Overlap Significance Tests ===")

    centered_path = RESULTS_DIR / "dirichlet_centered_full_rank.csv"
    uniform_path = RESULTS_DIR / "dirichlet_uniform_full_rank.csv"
    fixed_path = RESULTS_DIR / "top10_neural_tfs_prioritized.csv"

    methods = {}
    total_centered = total_uniform = total_fixed = 0

    if centered_path.exists():
        df_c = pd.read_csv(centered_path).drop_duplicates(subset="gene_id", keep="first")
        gene_col = "gene_id" if "gene_id" in df_c.columns else "gene_id_v6"
        score_col = "composite_score" if "composite_score" in df_c.columns \
            else "dirichlet_median_score"
        methods["centered"] = set(df_c.sort_values(score_col, ascending=False)[gene_col].head(10).values)
        total_centered = len(df_c)
    else:
        print(f"Warning: {centered_path} not found")
        methods["centered"] = set()

    if uniform_path.exists():
        df_u = pd.read_csv(uniform_path).drop_duplicates(subset="gene_id", keep="first")
        gene_col = "gene_id" if "gene_id" in df_u.columns else "gene_id_v6"
        score_col = "composite_score" if "composite_score" in df_u.columns \
            else "uniform_median_score"
        methods["uniform"] = set(df_u.sort_values(score_col, ascending=False)[gene_col].head(10).values)
        total_uniform = len(df_u)
    else:
        print(f"Warning: {uniform_path} not found")
        methods["uniform"] = set()

    if fixed_path.exists():
        # The PUBLISHED fixed-method shortlist (composite dual-track),
        # not a raw integrated_score top-10 from rank.csv.
        df_f = pd.read_csv(fixed_path)
        gene_col = "gene_id_v6" if "gene_id_v6" in df_f.columns else (
            "gene_id" if "gene_id" in df_f.columns else "gene_id_v6")
        methods["fixed"] = set(df_f[gene_col].astype(str))
        total_fixed = len(methods["fixed"])  # not used for N; N below
    else:
        print(f"Warning: {fixed_path} not found")
        methods["fixed"] = set()

    # Population size: the SHARED candidate universe - one row per gene in
    # rank.csv (never the row-exploded Dirichlet CSVs).
    rank_csv = RUN_DIR / "rank.csv"
    N = pd.read_csv(rank_csv)["gene_id"].nunique() if rank_csv.exists() else \
        max(total_centered, total_uniform)
    n = 10

    results = {"pairwise": {}, "three_way": {}, "binomial": {}, "overlaps": {}}

    method_names = list(methods.keys())
    print(f"\nMethods found: {method_names}")
    for m in method_names:
        print(f"  {m}: {len(methods[m])} top-10 genes")

    print("\n--- Pairwise Overlaps ---")
    for i in range(len(method_names)):
        for j in range(len(method_names)):
            if i == j:
                continue
            m1, m2 = method_names[i], method_names[j]
            overlap = methods[m1] & methods[m2]
            union = methods[m1] | methods[m2]
            k = len(overlap)
            jaccard = float(k / len(union)) if union else 0.0
            hg_p = hypergeometric_test(k, n, n, N)

            key = f"{m1}_vs_{m2}"
            results["pairwise"][key] = {
                "overlap_count": k,
                "jaccard": jaccard,
                "overlap_genes": sorted(list(overlap)),
                "hypergeometric_p": float(hg_p),
                "N_population": N,
            }
            if i < j:
                print(f"  {m1} vs {m2}: {k}/10 overlap (Jaccard={jaccard:.2f}), hypergeom p={hg_p:.4e}")

    for k, v in results["pairwise"].items():
        results["overlaps"][k] = {
            "count": v["overlap_count"],
            "jaccard": v["jaccard"],
            "p_value": v["hypergeometric_p"],
        }

    print("\n--- Three-way Overlap ---")
    if len(methods) == 3:
        three_way = methods.get("centered", set()) & methods.get("uniform", set()) & methods.get("fixed", set())
        k3 = len(three_way)
        # Proper 3-set null: P(gene in all 3 random n-subsets) = (n/N)^3;
        # count ~ Binomial(N, (n/N)^3). The previous pairwise-null call
        # (hypergeometric with K=n) was conservative by ~8 orders of
        # magnitude for small k3.
        p3 = three_way_null_p(k3, n, N)
        results["three_way"] = {
            "overlap_count": k3,
            "overlap_genes": sorted(list(three_way)),
            "null_p_per_gene": float((n / N) ** 3),
            "binomial_three_set_p": float(p3),
            "N_population": N,
        }
        results["overlaps"]["three_way"] = {
            "count": k3,
            "p_value": float(p3),
        }
        print(f"  Three-way overlap: {k3}/10, 3-set binomial p={p3:.4e}")

    print("\n--- Binomial Test ---")
    unique_pairs = [(method_names[i], method_names[j]) for i in range(len(method_names)) for j in range(i+1, len(method_names))]
    total_possible_pairs = len(unique_pairs)
    total_overlap_count = sum(
        results["pairwise"][f"{m1}_vs_{m2}"]["overlap_count"]
        for m1, m2 in unique_pairs
        if f"{m1}_vs_{m2}" in results["pairwise"]
    )
    max_possible = total_possible_pairs * n
    if max_possible > 0:
        binom_p = _binom_p(total_overlap_count, max_possible, 1.0 / N * n)
        results["binomial"] = {
            "total_overlaps": total_overlap_count,
            "max_possible": max_possible,
            "observed_rate": total_overlap_count / max_possible if max_possible > 0 else 0,
            "expected_rate_under_null": n / N,
            "binomial_p_greater": float(binom_p),
        }
        print(f"  Total overlaps: {total_overlap_count}/{max_possible}, binom p={binom_p:.4e}")
    else:
        results["binomial"] = {"error": "no pairs"}

    out_path = RESULTS_DIR / "overlap_significance.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
