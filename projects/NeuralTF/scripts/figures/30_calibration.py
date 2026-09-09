"""Empirical decile calibration diagram for the integrated evidence score.

Single-panel evaluation of rank discrimination and score calibration:
- Validated neural TF rate across score deciles (D1 lowest to D10 highest)
- Wilson 95% binomial confidence intervals
- Genome-wide background prevalence reference line
- Top-decile enrichment fold-change annotation
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json


def build():
    data_path = RES / "calibration_stats.json"
    if not data_path.exists():
        raise FileNotFoundError(f"{data_path} missing — run scripts/stats/calibration.py first")

    with open(data_path) as f:
        data = json.load(f)

    stats_list = data.get("bin_stats", [])
    if not stats_list:
        raise ValueError("calibration_stats.json carries no bin_stats")
    bin_df = pd.DataFrame(stats_list)

    # Order deciles from D1 (lowest) to D10 (highest score)
    bin_df = bin_df.sort_values("mean_score", ascending=True).reset_index(drop=True)
    decile_labels = [f"D{i+1}\n({bin_df.iloc[i]['mean_score']:.2f})" for i in range(len(bin_df))]
    x = np.arange(len(bin_df))
    observed = bin_df["empirical_positive_rate"].to_numpy(dtype=float)
    counts = bin_df["n_candidates"].to_numpy(dtype=float)
    prevalence = float(data["prevalence"])

    fig, ax = plt.subplots(figsize=(W_15COL, 3.8), dpi=500)

    # Bars: Highlight bins enriched above background prevalence
    colors = [C_A if obs >= prevalence else "#D0D7DE" for obs in observed]
    bars = ax.bar(x, observed * 100, color=colors, width=0.62, edgecolor="none", zorder=3)

    # Wilson 95% confidence intervals per decile
    z = 1.96
    p_hat = np.clip(observed, 0, 1)
    denom = 1 + z**2 / counts
    half = z * np.sqrt(p_hat * (1 - p_hat) / counts + z**2 / (4 * counts**2)) / denom
    ax.errorbar(x, observed * 100, yerr=half * 100, fmt="none", ecolor="#333333",
                elinewidth=0.9, capsize=3, zorder=4, label="Wilson 95% CI")

    # Prevalence reference line
    ax.axhline(y=prevalence * 100, color=C_HL, lw=1.1, linestyle="--", zorder=2,
                label=f"Genome prevalence ({prevalence*100:.2f}%)")

    # Top decile callout annotation
    top_rate = observed[-1] * 100
    enrichment = observed[-1] / prevalence if prevalence > 0 else 0
    ax.annotate(f"{top_rate:.1f}%\n({enrichment:.1f}\u00d7 enrichment)",
                xy=(x[-1], top_rate + half[-1] * 100),
                xytext=(x[-1] - 1.2, top_rate + half[-1] * 100 + 0.6),
                arrowprops=dict(arrowstyle="->", color="#333333", lw=0.7),
                fontsize=7.2, fontweight="bold", ha="center", color="#222222")

    ax.set_xticks(x)
    ax.set_xticklabels(decile_labels, fontsize=6.8)
    ax.set_xlabel("Integrated Score Decile (Mean Score)", fontsize=8, fontweight="bold")
    ax.set_ylabel("RNAi Validation Rate (%)", fontsize=8, fontweight="bold")
    ax.set_title("Rank Discrimination and Score Calibration Across Transcriptome Deciles",
                 fontweight="bold", fontsize=8.5, pad=8)
    ax.legend(loc="upper left", frameon=False, fontsize=7)
    ax.set_ylim(-0.2, max(observed * 100 + half * 100) * 1.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.subplots_adjust(left=0.14, right=0.96, top=0.90, bottom=0.15)
    save(fig, "30_calibration")


if __name__ == "__main__":
    build()
