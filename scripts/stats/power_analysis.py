#!/usr/bin/env python
"""Power analysis and convergence diagnostics.

Analyzes convergence of Dirichlet weight draws and statistical power
of the permutation test framework.

Usage:
    python scripts/stats/power_analysis.py --n-draws 500 --n-perm-test 200 --seed 42
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


def integrated_score(S, W):
    """Compute integrated score with missing-data renormalization."""
    mask = ~np.isnan(S)
    if not mask.any():
        return 0.0
    S_filled = np.where(np.isnan(S), 0.0, S)
    w_masked = W[mask]
    return np.sum(S_filled[mask] * w_masked) / w_masked.sum()


def simulate_scores(n_candidates=249, n_streams=8, rng=None):
    """Simulate evidence stream scores for power analysis."""
    if rng is None:
        rng = np.random.default_rng(42)

    scores = rng.uniform(0, 1, size=(n_candidates, n_streams))

    n_true = 20
    for i in range(n_true):
        scores[i, :] = rng.uniform(0.5, 1.0, size=n_streams)

    missing_rate = 0.15
    missing_mask = rng.random(size=scores.shape) < missing_rate
    scores[missing_mask] = np.nan

    return scores, n_true


def convergence_analysis(scores, weights, n_draws_list, rng):
    """Analyze convergence of integrated scores as function of Dirichlet draws."""
    n_candidates = scores.shape[0]
    results = {}

    for n_draws in n_draws_list:
        rank_stabilities = []
        score_stabilities = []

        for _ in range(20):
            draw_weights = rng.dirichlet(np.ones(len(weights)), size=n_draws)
            mean_weights = draw_weights.mean(axis=0)
            mean_weights = mean_weights / mean_weights.sum()

            rep_scores = np.array([integrated_score(row, mean_weights) for row in scores])
            rank_stabilities.append(np.corrcoef(
                np.argsort(-np.array([integrated_score(row, weights) for row in scores])),
                np.argsort(-rep_scores)
            )[0, 1])
            score_stabilities.append(np.corrcoef(
                np.array([integrated_score(row, weights) for row in scores]),
                rep_scores
            )[0, 1])

        results[n_draws] = {
            "rank_stability_mean": np.mean(rank_stabilities),
            "rank_stability_std": np.std(rank_stabilities),
            "score_stability_mean": np.mean(score_stabilities),
            "score_stability_std": np.std(score_stabilities),
        }

    return results


def power_permutation_test(scores, weights, n_true, n_perms_list, rng):
    """Estimate power of permutation test for different numbers of permutations."""
    true_integrated = np.array([integrated_score(row, weights) for row in scores])

    power_results = {}
    for n_perms in n_perms_list:
        sig_count = 0
        for _ in range(50):
            null_scores = []
            for _ in range(n_perms):
                perm_idx = rng.permutation(scores.shape[0])
                perm_scores = np.array([integrated_score(row, weights) for row in scores[perm_idx]])
                null_scores.append(np.max(perm_scores))
            null_dist = np.array(null_scores)
            p_empirical = (np.sum(null_dist >= true_integrated[0]) + 1) / (n_perms + 1)
            if p_empirical < 0.05:
                sig_count += 1

        power_results[n_perms] = {
            "estimated_power": sig_count / 50,
            "n_perms": n_perms,
        }

    return power_results


def main():
    parser = argparse.ArgumentParser(description="Power analysis and convergence diagnostics")
    parser.add_argument("--n-draws", type=int, default=500, help="Max Dirichlet draws to test")
    parser.add_argument("--n-perm-test", type=int, default=200, help="Max permutations for power test")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    print("=== Power Analysis & Convergence Diagnostics ===")

    rng = np.random.default_rng(args.seed)
    scores, n_true = simulate_scores(rng=rng)
    weights = np.array([0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2])
    weights = weights / weights.sum()

    print(f"Simulated data: {scores.shape[0]} candidates, {scores.shape[1]} streams")
    print(f"True signals: {n_true}")

    print("\n--- Convergence Analysis ---")
    draws_list = [10, 25, 50, 100, 200, 500]
    draws_list = [d for d in draws_list if d <= args.n_draws]
    conv_results = convergence_analysis(scores, weights, draws_list, rng)

    for n_draws, res in conv_results.items():
        print(f"  Draws={n_draws:>4d}: rank_stability={res['rank_stability_mean']:.4f} "
              f"± {res['rank_stability_std']:.4f}, "
              f"score_stability={res['score_stability_mean']:.4f}")

    print("\n--- Permutation Test Power ---")
    perms_list = [50, 100, 200]
    perms_list = [p for p in perms_list if p <= args.n_perm_test]
    power_results = power_permutation_test(scores, weights, n_true, perms_list, rng)

    for n_perms, res in power_results.items():
        print(f"  Perms={n_perms:>4d}: estimated_power={res['estimated_power']:.2f}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    draws_vals = sorted(conv_results.keys())
    rank_means = [conv_results[d]["rank_stability_mean"] for d in draws_vals]
    rank_stds = [conv_results[d]["rank_stability_std"] for d in draws_vals]
    ax.errorbar(draws_vals, rank_means, yerr=rank_stds, marker="o", linewidth=2, capsize=4)
    ax.set_xlabel("Number of Dirichlet Draws")
    ax.set_ylabel("Spearman Rank Stability")
    ax.set_title("Convergence of Dirichlet Weight Averaging")
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0.8, 1.02])

    ax = axes[1]
    perms_vals = sorted(power_results.keys())
    power_vals = [power_results[p]["estimated_power"] for p in perms_vals]
    ax.plot(perms_vals, power_vals, "ro-", linewidth=2, markersize=8)
    ax.set_xlabel("Number of Permutations")
    ax.set_ylabel("Estimated Power")
    ax.set_title("Permutation Test Power vs. Permutations")
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.1])
    ax.axhline(y=0.8, color="gray", linestyle="--", alpha=0.5, label="80% power")
    ax.legend()

    plt.tight_layout()
    fig_path = FIG_DIR / "convergence_analysis.png"
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)
    print(f"\nSaved: {fig_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
