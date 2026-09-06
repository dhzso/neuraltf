"""Calibration — decile rank-discrimination diagram for the integrated score.

The integrated score is NOT a probability, so this is a rank-
discrimination plot (empirical RNAi-validated fraction per score decile,
monotonicity = discrimination power), NOT a probability reliability
diagram. The misleading perfect-calibration diagonal and the
score-vs-fraction pairing are removed (the score is not P(positive)).
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json

def build():
    data_path = RES / "calibration_stats.json"
    if not data_path.exists():
        raise FileNotFoundError(f"{data_path} missing - run scripts/stats/calibration.py first")

    with open(data_path) as f:
        data = json.load(f)

    fig, ax = plt.subplots(figsize=(6.5, 5))

    # Decile bins in score space: empirical positive rate per decile
    stats_list = data.get("bin_stats", [])
    if not stats_list:
        raise ValueError("calibration_stats.json carries no bin_stats")
    bin_df = pd.DataFrame(stats_list)
    # bin index 0 = lowest decile .. n-1 = highest (calibration.py orders
    # deciles ascending then inverts label 9=highest; re-derive here from
    # mean_score ordering)
    bin_df = bin_df.sort_values("mean_score", ascending=True).reset_index(drop=True)
    x = np.arange(len(bin_df))
    observed = bin_df["empirical_positive_rate"].to_numpy(dtype=float)
    counts = bin_df["n_candidates"].to_numpy(dtype=float)
    prevalence = float(data["prevalence"])

    bars = ax.bar(x, observed, color=C_A, alpha=0.85, width=0.7,
                  label="Empirical RNAi-validated fraction")
    ax.axhline(y=prevalence, color=C_HL, lw=1.2, linestyle="--",
               label=f"Cohort prevalence ({prevalence:.4f})")

    # Wilson 95% CI per bin (prevalence is tiny; Wilson stays in [0,1])
    z = 1.96
    p_hat = np.clip(observed, 0, 1)
    denom = 1 + z**2 / counts
    center = (p_hat + z**2 / (2 * counts)) / denom
    half = z * np.sqrt(p_hat * (1 - p_hat) / counts + z**2 / (4 * counts**2)) / denom
    ax.errorbar(x, observed, yerr=half, fmt="none", ecolor="#333333",
                elinewidth=1, capsize=3)

    # secondary axis: bin counts
    ax2 = ax.twinx()
    ax2.plot(x, counts, color=C_NEURAL, lw=1.5, marker=".", markersize=5,
             alpha=0.7, label="Candidates per decile")
    ax2.set_ylabel("Candidates per decile", color="#999999")
    ax2.tick_params(axis="y", labelcolor="#999999")

    ax.set_xlabel("Integrated-score decile (low → high)")
    ax.set_ylabel("Fraction RNAi-validated in decile")
    ax.set_title("Rank discrimination by score decile\n"
                 "Monotone rise = the score enriches validated neural TFs "
                 "(not a probability calibration)",
                 fontweight="bold", pad=8)
    ax.legend(loc="upper left", fontsize=8)
    ax.spines["top"].set_visible(False)

    fig.tight_layout()
    save(fig, "30_calibration")

if __name__ == "__main__":
    build()
