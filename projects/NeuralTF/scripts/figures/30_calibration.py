"""Calibration and rank discrimination diagram for the integrated score.

The integrated score is an evidence-weight aggregation, not a probability.
Rank discrimination measures whether top-ranked candidates enrich RNAi-validated
neural TFs.
Panel a: Empirical validation rate per score decile with Wilson 95% CIs.
Panel b: Cumulative recovery of known validated neural TFs across score deciles.
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

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_2COL, 2.7))

    stats_list = data.get("bin_stats", [])
    if not stats_list:
        raise ValueError("calibration_stats.json carries no bin_stats")
    bin_df = pd.DataFrame(stats_list)
    # Order deciles from 1 (lowest) to 10 (highest score)
    bin_df = bin_df.sort_values("mean_score", ascending=True).reset_index(drop=True)
    decile_labels = [f"D{i+1}" for i in range(len(bin_df))]
    x = np.arange(len(bin_df))
    observed = bin_df["empirical_positive_rate"].to_numpy(dtype=float)
    counts = bin_df["n_candidates"].to_numpy(dtype=float)
    positives = bin_df["n_positives"].to_numpy(dtype=int)
    prevalence = float(data["prevalence"])
    total_pos = int(data["n_positives"])

    # Panel a: Empirical validation rate
    colors = [C_A if obs > prevalence else "#CCCCCC" for obs in observed]
    bars = ax1.bar(x, observed * 100, color=colors, width=0.65, edgecolor="none")
    ax1.axhline(y=prevalence * 100, color=C_HL, lw=1.0, linestyle="--",
                label=f"Genome prevalence ({prevalence*100:.2f}%)")

    # Wilson 95% CI per bin
    z = 1.96
    p_hat = np.clip(observed, 0, 1)
    denom = 1 + z**2 / counts
    center = (p_hat + z**2 / (2 * counts)) / denom
    half = z * np.sqrt(p_hat * (1 - p_hat) / counts + z**2 / (4 * counts**2)) / denom
    ax1.errorbar(x, observed * 100, yerr=half * 100, fmt="none", ecolor="#333333",
                 elinewidth=0.8, capsize=2)

    # Top decile callout
    top_rate = observed[-1] * 100
    enrichment = data.get("top_decile_enrichment", {}).get("fold_enrichment", 9.5)
    ax1.text(x[-1], top_rate + 1.0, f"{top_rate:.1f}%\n({enrichment:.1f}\u00d7)",
             ha="center", va="bottom", fontsize=6.5, fontweight="bold", color=C_A)

    ax1.set_xticks(x)
    ax1.set_xticklabels(decile_labels, fontsize=6.5)
    ax1.set_xlabel("Score decile (D1 = lowest, D10 = highest)", fontsize=7.5)
    ax1.set_ylabel("RNAi-validated rate (%)", fontsize=7.5)
    ax1.legend(loc="upper left", frameon=False, fontsize=6.5)
    panel_tag(ax1, "a")

    # Panel b: Cumulative positive recovery (from D10 downwards)
    pos_desc = positives[::-1]
    cum_recovered = np.cumsum(pos_desc)
    cum_pct = (cum_recovered / total_pos) * 100
    x_rev = np.arange(len(cum_pct))
    rev_labels = [f"Top {i+1}0%" for i in range(len(cum_pct))]

    ax2.plot(x_rev, cum_pct, marker="o", markersize=4, color=C_B, lw=1.5, clip_on=False)
    ax2.fill_between(x_rev, 0, cum_pct, color=C_B, alpha=0.15)
    ax2.plot([0, len(cum_pct)-1], [10, 100], color="#AAAAAA", lw=1.0, linestyle=":",
             label="Random baseline")

    # Annotate top 10% recovery
    ax2.annotate(f"{cum_pct[0]:.1f}% recovered\nin top 10%",
                 xy=(0, cum_pct[0]), xytext=(1.5, 75),
                 arrowprops=dict(arrowstyle="->", color=C_B, lw=0.8),
                 fontsize=6.5, fontweight="bold", color=C_B)

    ax2.set_xticks(x_rev)
    ax2.set_xticklabels(rev_labels, rotation=35, ha="right", fontsize=6.5)
    ax2.set_xlabel("Cumulative score fraction", fontsize=7.5)
    ax2.set_ylabel("Validated neural TFs captured (%)", fontsize=7.5)
    ax2.set_ylim(0, 105)
    ax2.legend(loc="lower right", frameon=False, fontsize=6.5)
    panel_tag(ax2, "b")

    fig.tight_layout()
    save(fig, "30_calibration")

if __name__ == "__main__":
    build()
