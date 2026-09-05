#!/usr/bin/env python
"""Dirichlet draw convergence + permutation-resolution diagnostics.

2026-09-04 statistical redesign (audit findings):

Panel 1 (kept, metric fixed): how many Dirichlet draws are needed before
the rank list stabilizes? The previous code Spearman-correlated two
argsort ORDER VECTORS (arrays of gene indices in rank order), which is
NOT a rank correlation of the rankings - it returns ~0 for any pair of
non-identical orderings (and ~0 even for perfectly reversed rankings).
The correct metric: convert each ordering to per-gene RANK vectors and
Spearman-correlate those (equivalently, Pearson on the ranks; identical
to Spearman of the two score vectors).

Panel 2 (replaced): the previous "power" panel was tautological - it
defined power as P(the top-scoring gene reaches p<0.05) under a label
permutation where p could only take two values, both < 0.05 for n>=40;
observed power was 1.0/1.0/1.0 by construction. It is replaced by a
permutation-RESOLUTION analysis: the empirical p-value granularity of
the production permutation test as a function of n_perm -
p_floor = 1/(n_perm+1) - and the minimum detectable p, which is the
honest way to choose n_perm for the real test.

Writes the inputs consumed by figure 29:
  results/convergence_draws.csv
  results/permutation_resolution.csv

Usage:
    python scripts/stats/power_analysis.py
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DRAW_PATH = RESULTS_DIR / "dirichlet_centered_draw_scores.csv"


def load_draw_matrix() -> tuple[np.ndarray, np.ndarray] | None:
    if not DRAW_PATH.exists():
        return None
    draws = pd.read_csv(DRAW_PATH).drop_duplicates(subset="gene_id", keep="first")
    gene_col = "gene_id" if "gene_id" in draws.columns else draws.columns[0]
    draw_cols = [c for c in draws.columns if c.startswith("draw_")]
    if len(draw_cols) < 50:
        return None
    return draws[gene_col].astype(str).to_numpy(), \
        draws[draw_cols].to_numpy(dtype=float)


def _order_to_ranks(order: np.ndarray) -> np.ndarray:
    """Convert an ordering (positions -> gene index) into per-gene ranks
    (gene index -> rank). argsort outputs positions; invert them."""
    ranks = np.empty(len(order), dtype=float)
    ranks[order] = np.arange(1, len(order) + 1)
    return ranks


def convergence_from_draws(mat: np.ndarray, sizes, rng, n_rep=20):
    """Spearman stability of prefix-median RANK vectors vs the full-draw
    rank vector (correct metric: ranks correlated, not order vectors)."""
    full_ranks = _order_to_ranks(np.argsort(-np.median(mat, axis=1), kind="stable"))
    results = {}
    for size in sizes:
        if size > mat.shape[1]:
            continue
        rhos = []
        for _ in range(n_rep):
            take = rng.choice(mat.shape[1], size=size, replace=False)
            prefix_ranks = _order_to_ranks(
                np.argsort(-np.median(mat[:, take], axis=1), kind="stable"))
            rho = stats.spearmanr(full_ranks, prefix_ranks).statistic
            rhos.append(0.0 if np.isnan(rho) else float(rho))
        results[size] = {
            "n_draws": size,
            "spearman_vs_full": float(np.mean(rhos)),
            "spearman_std": float(np.std(rhos)),
            "n_rep": n_rep,
        }
    return results


def permutation_resolution(perms_list):
    """Empirical p-value resolution of the (b+1)/(n+1) permutation test
    as a function of n_perm: the smallest achievable p (all-null-lost
    case) and the smallest p that can reach the 0.05 threshold.

    Replaces the tautological 'power' panel: under the exchangeable null,
    p_floor = 1/(n+1) is the minimum possible p-value, so n_perm must
    satisfy 1/(n_perm+1) < alpha for the test to have ANY capacity to
    reject at level alpha - the honest n_perm guidance.
    """
    results = {}
    for n_perm in perms_list:
        p_floor = 1.0 / (n_perm + 1)
        results[n_perm] = {
            "n_perm": n_perm,
            "min_detectable_p": p_floor,
            "resolves_p05": bool(p_floor < 0.05),
            "resolves_p01": bool(p_floor < 0.01),
        }
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Draw convergence + permutation resolution diagnostics"
    )
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    print("=== Draw Convergence + Permutation Resolution (real data) ===")
    rng = np.random.default_rng(args.seed)

    loaded = load_draw_matrix()
    if loaded is None:
        print(
            f"ERROR: {DRAW_PATH} missing or too small.\n"
            "Run projects/NeuralTF/scripts/dirichlet_centered.py first - "
            "it persists the per-candidate draw-score matrix this analysis "
            "consumes. (The previous version plotted simulated data.)"
        )
        return 1
    genes, mat = loaded
    print(f"Draw matrix: {mat.shape[0]} candidates x {mat.shape[1]} draws")

    # --- Panel 1: convergence from the real draws (correct metric) ------
    sizes = [10, 25, 50, 100, 250, 500, mat.shape[1]]
    conv = convergence_from_draws(mat, sizes, rng)
    conv_df = pd.DataFrame([conv[s] for s in sorted(conv)])
    conv_path = RESULTS_DIR / "convergence_draws.csv"
    conv_df.to_csv(conv_path, index=False)
    print(f"\nConvergence (Spearman of prefix RANKS vs full {mat.shape[1]} draws):")
    for _, r in conv_df.iterrows():
        print(f"  draws={int(r['n_draws']):>4}: rho={r['spearman_vs_full']:.4f} "
              f"+- {r['spearman_std']:.4f}")
    print(f"Saved: {conv_path}")

    # --- Panel 2: permutation resolution guidance ------------------------
    perms_list = [10, 20, 30, 50, 100, 200, 500, 1000]
    res = permutation_resolution(perms_list)
    res_df = pd.DataFrame([res[p] for p in sorted(res)])
    res_path = RESULTS_DIR / "permutation_resolution.csv"
    res_df.to_csv(res_path, index=False)
    print(f"\nPermutation resolution (add-one p-value granularity):")
    for _, r in res_df.iterrows():
        print(f"  n_perm={int(r['n_perm']):>4}: min p={r['min_detectable_p']:.5f} "
              f"(reject at 0.05: {r['resolves_p05']})")
    print(f"Saved: {res_path}")

    # Retire the tautological panel's file so no consumer reads stale
    # "power 1.0" values.
    stale = RESULTS_DIR / "permutation_power.csv"
    if stale.exists():
        stale.unlink()
        print(f"Removed obsolete tautological panel file: {stale.name}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
