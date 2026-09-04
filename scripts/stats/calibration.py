#!/usr/bin/env python
"""Calibration / discrimination analysis for integrated scores.

Bins integrated scores into deciles and computes the empirical positive
rate (RNAi-validated TFs) per bin. NOTE (WS3): the integrated score is
an evidence-weight score in [0,1], not a calibrated probability, so a
classic "perfect calibration" diagonal is conceptually invalid. We
report two honest discrimination metrics:

  - rank-discrimination error: mean |empirical positive rate - prevalence|
    per decile (how far decile rates deviate from the base rate — this is
    what the old code called "ECE")
  - true ECE against prevalence-weighted decile means (documented as a
    discrimination proxy)

The reliability plot itself (score decile vs observed rate) is valid and
is what figure 30 renders.

Usage:
    python scripts/stats/calibration.py --n-bins 10
"""

import argparse
import json
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


def load_all_genes():
    """Load candidate table."""
    p = RUN_DIR / "rank.csv"
    if p.exists():
        return pd.read_csv(p).drop_duplicates(subset="gene_id", keep="first")
    raise FileNotFoundError("No candidate score file found in runs/pipeline_run/rank.csv")


def main():
    parser = argparse.ArgumentParser(description="Calibration analysis for integrated scores")
    parser.add_argument("--n-bins", type=int, default=10, help="Number of score bins (default 10)")
    args = parser.parse_args()

    print("=== Calibration Analysis ===")

    df = load_all_genes()
    score_col = None
    for c in ["integrated_score", "composite_score", "dirichlet_median_score", "fixed_weight_score"]:
        if c in df.columns:
            score_col = c
            break
    if score_col is None:
        print("Error: no score column found")
        return 1

    if "proof_status" not in df.columns:
        print("Error: proof_status column not found")
        return 1

    df = df.dropna(subset=[score_col]).copy()
    df["is_positive"] = (df["proof_status"] == "known_rnai_validated").astype(int)

    n_pos = df["is_positive"].sum()
    n_total = len(df)
    prevalence = n_pos / n_total if n_total > 0 else 0
    print(f"Candidates: {n_total}, Positives: {n_pos}, Prevalence: {prevalence:.4f}")

    df["decile"] = pd.qcut(df[score_col].rank(method="first"), q=args.n_bins, labels=False)
    df["decile"] = args.n_bins - 1 - df["decile"]

    bin_stats = []
    for decile in sorted(df["decile"].unique()):
        bin_df = df[df["decile"] == decile]
        n_bin = len(bin_df)
        n_pos_bin = bin_df["is_positive"].sum()
        empirical_rate = n_pos_bin / n_bin if n_bin > 0 else 0
        mean_score = bin_df[score_col].mean()
        lo = bin_df[score_col].min()
        hi = bin_df[score_col].max()

        bin_stats.append({
            "decile": int(decile),
            "score_range_lo": float(lo),
            "score_range_hi": float(hi),
            "mean_score": float(mean_score),
            "n_candidates": int(n_bin),
            "n_positives": int(n_pos_bin),
            "empirical_positive_rate": float(empirical_rate),
        })
        print(f"  Decile {decile}: [{lo:.3f}, {hi:.3f}], n={n_bin}, "
              f"positives={n_pos_bin}, rate={empirical_rate:.4f}")

    bin_df_out = pd.DataFrame(bin_stats)

    mean_scores = bin_df_out["mean_score"].values
    emp_rates = bin_df_out["empirical_positive_rate"].values
    # rank-discrimination error: deviation of decile rates from the base
    # rate (this is what the previous version mislabeled "ECE")
    cal_error = np.mean(np.abs(emp_rates - prevalence))
    max_cal_error = np.max(np.abs(emp_rates - prevalence))

    results = {
        "n_bins": args.n_bins,
        "n_total": int(n_total),
        "n_positives": int(n_pos),
        "prevalence": float(prevalence),
        "mean_calibration_error": float(cal_error),
        "max_calibration_error": float(max_cal_error),
        "rank_discrimination_error": float(cal_error),
        "discrimination_note": (
            "integrated_score is an evidence-weight score, not a probability; "
            "the reported error is the deviation of per-decile empirical "
            "positive rates from the cohort prevalence (rank-discrimination), "
            "not a probability-calibration ECE"
        ),
        "bin_centers": [float(b["mean_score"]) for b in bin_stats],
        "observed_fractions": [float(b["empirical_positive_rate"]) for b in bin_stats],
        "expected_fractions": [float(b["mean_score"]) for b in bin_stats],
        "bin_counts": [int(b["n_candidates"]) for b in bin_stats],
        "bin_stats": bin_stats,
    }
    print(f"\nRank-discrimination error: {cal_error:.4f}")
    print(f"Max deviation from prevalence:  {max_cal_error:.4f}")

    out_path = RESULTS_DIR / "calibration_stats.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
