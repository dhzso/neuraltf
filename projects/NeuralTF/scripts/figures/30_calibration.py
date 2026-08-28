"""Calibration — reliability diagram for integrated score bins."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import json

def build():
    data_path = RES / "calibration_stats.json"
    with open(data_path) as f:
        data = json.load(f)

    fig, ax = plt.subplots(figsize=(6, 5))

    bins = np.array(data["bin_centers"])
    observed = np.array(data["observed_fractions"])
    expected = np.array(data["expected_fractions"])
    counts = np.array(data["bin_counts"])

    ax.plot([0, 1], [0, 1], color="#999999", lw=1, linestyle="--", label="Perfect calibration")
    ax.plot(bins, observed, color=C_A, lw=2, marker="o", markersize=6, label="Observed fraction")

    # Error bars for binomial proportion
    se = np.sqrt(observed * (1 - observed) / counts)
    ax.errorbar(bins, observed, yerr=1.96 * se, fmt="none", ecolor=C_A, elinewidth=1, capsize=3)

    # Bar width for bin counts (secondary axis)
    ax2 = ax.twinx()
    bar_width = 0.04
    ax2.bar(bins, counts, width=bar_width, color=C_NEURAL, alpha=0.3, label="Bin count")
    ax2.set_ylabel("Number of candidates per bin", color="#999999")
    ax2.tick_params(axis="y", labelcolor="#999999")

    # ECE
    ece = data.get("ece", None)
    if ece is not None:
        ax.text(0.05, 0.95, f"ECE = {ece:.3f}", transform=ax.transAxes,
                fontsize=9, fontweight="bold", va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=C_HL, alpha=0.8))

    ax.set_xlabel("Predicted score bin (expected fraction of true positives)")
    ax.set_ylabel("Observed fraction of true positives")
    ax.set_title("Reliability diagram: integrated score calibration\n"
                 "Closer to diagonal = better calibrated",
                 fontweight="bold", pad=8)
    ax.legend(loc="upper left", fontsize=8)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.spines["top"].set_visible(False)

    fig.tight_layout()
    save(fig, "30_calibration")

if __name__ == "__main__":
    build()
