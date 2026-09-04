#!/usr/bin/env python
"""Power analysis and convergence diagnostics (REAL draw matrices).

Panel 1 — Convergence: how many Dirichlet draws are needed before the
rank list stabilizes? Computed from the pipeline's actual persisted draw
matrices (results/dirichlet_centered_draw_scores.csv): for prefix sizes
10..1000, Spearman-correlate the prefix-median ranking against the
full-1000-draw ranking (20 bootstrap repetitions of prefix choice).

Panel 2 — Power: Monte-Carlo estimate of the permutation test's power as
a function of n_perm, using the real observed score distribution and the
same (count+1)/(n+1) p-value estimator the permutation test uses.

Writes the inputs consumed by figure 29:
  results/convergence_draws.csv
  results/permutation_power.csv
(and no standalone figure — the numbered publication script renders it).

Usage:
    python scripts/stats/power_analysis.py [--n-perm-test 200] [--seed 42]
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


def convergence_from_draws(mat: np.ndarray, sizes, rng, n_rep=20):
    """Spearman stability of prefix-median rankings vs the full-draw ranking."""
    full_order = np.argsort(-np.median(mat, axis=1))
    results = {}
    for size in sizes:
        if size > mat.shape[1]:
            continue
        rhos = []
        for _ in range(n_rep):
            take = rng.choice(mat.shape[1], size=size, replace=False)
            prefix_order = np.argsort(-np.median(mat[:, take], axis=1))
            rho = stats.spearmanr(full_order, prefix_order).statistic
            rhos.append(0.0 if np.isnan(rho) else float(rho))
        results[size] = {
            "n_draws": size,
            "spearman_vs_full": float(np.mean(rhos)),
            "spearman_std": float(np.std(rhos)),
            "n_rep": n_rep,
        }
    return results


def power_vs_nperm(rank: pd.DataFrame, score_col: str, perms_list, rng,
                   n_rep=30):
    """MC power of the permutation test at different n_perm.

    Null: score labels randomly permuted across genes (assignment null);
    p = (count+1)/(n+1) exactly as in permutation_test_full.py. Power is
    the fraction of repetitions in which the top-scoring gene reaches
    p < 0.05 — i.e., the test's ability to flag a genuine top candidate.
    """
    scores = pd.to_numeric(rank[score_col], errors="coerce").dropna().to_numpy()
    results = {}
    for n_perm in perms_list:
        hits = 0
        for _ in range(n_rep):
            nulls = rng.permutation(scores)[:n_perm]
            # one-sided: how often a random assignment reaches the top score
            p_emp = (np.sum(nulls >= scores.max()) + 1) / (n_perm + 1)
            if p_emp < 0.05:
                hits += 1
        results[n_perm] = {"n_perm": n_perm, "power": hits / n_rep}
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Power analysis + convergence diagnostics on real draws"
    )
    parser.add_argument("--n-perm-test", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    print("=== Power Analysis & Convergence Diagnostics (real data) ===")
    rng = np.random.default_rng(args.seed)

    loaded = load_draw_matrix()
    if loaded is None:
        print(
            f"ERROR: {DRAW_PATH} missing or too small.\n"
            "Run projects/NeuralTF/scripts/dirichlet_centered.py first — "
            "it persists the per-candidate draw-score matrix this analysis "
            "consumes. (The previous version plotted simulated data.)"
        )
        return 1
    genes, mat = loaded
    print(f"Draw matrix: {mat.shape[0]} candidates x {mat.shape[1]} draws")

    # --- Panel 1: convergence from the real draws -------------------------
    sizes = [10, 25, 50, 100, 250, 500, mat.shape[1]]
    conv = convergence_from_draws(mat, sizes, rng)
    conv_df = pd.DataFrame([conv[s] for s in sorted(conv)])
    conv_path = RESULTS_DIR / "convergence_draws.csv"
    conv_df.to_csv(conv_path, index=False)
    print(f"\nConvergence (Spearman of prefix ranking vs full {mat.shape[1]} draws):")
    for _, r in conv_df.iterrows():
        print(f"  draws={r['n_draws']:>4}: rho={r['spearman_vs_full']:.4f} "
              f"± {r['spearman_std']:.4f}")
    print(f"Saved: {conv_path}")

    # --- Panel 2: permutation power on the real score vector ---------------
    rank_path = RUN_DIR / "rank.csv"
    if not rank_path.exists():
        print(f"\n[skip power panel] {rank_path} missing (run the pipeline first)")
        return 0
    rank = pd.read_csv(rank_path).drop_duplicates(subset="gene_id", keep="first")
    score_col = "integrated_score" if "integrated_score" in rank.columns \
        else rank.columns[-1]
    perms_list = [50, 100, 200]
    perms_list = [p for p in perms_list if p <= args.n_perm_test]
    power = power_vs_nperm(rank, score_col, perms_list, rng)
    power_df = pd.DataFrame([power[p] for p in sorted(power)])
    power_path = RESULTS_DIR / "permutation_power.csv"
    power_df.to_csv(power_path, index=False)
    print(f"\nPermutation power (top candidate reaches p<0.05):")
    for _, r in power_df.iterrows():
        print(f"  n_perm={r['n_perm']:>4}: power={r['power']:.2f}")
    print(f"Saved: {power_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
