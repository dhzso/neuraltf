"""Dirichlet draw convergence analysis.

Single-panel plot demonstrating convergence of candidate ranking stability as a function
of Dirichlet draws (k = 40) against the full 1,000-draw reference ranking.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def build():
    conv_path = RES / "convergence_draws.csv"
    if not conv_path.exists():
        raise FileNotFoundError(
            f"{conv_path} missing — run scripts/stats/power_analysis.py first"
        )
    conv = pd.read_csv(conv_path)

    fig, ax = plt.subplots(figsize=(W_15COL, 3.6), dpi=500)

    x = conv["n_draws"].values
    y = conv["spearman_vs_full"].values
    err = conv["spearman_std"].values if "spearman_std" in conv.columns else None

    if err is not None:
        ax.errorbar(x, y, yerr=err, color=C_A, lw=1.6, marker="o",
                    markersize=4.5, capsize=3, ecolor="#555555",
                    label="Observed rank correlation (mean \u00b1 s.d.)")
    else:
        ax.plot(x, y, color=C_A, lw=1.6, marker="o", markersize=4.5,
                label="Observed rank correlation (Spearman $r_s$)")

    ax.axhline(y=0.99, color="#4A7C59", lw=0.9, linestyle="--", alpha=0.85,
               label="High fidelity threshold ($r_s = 0.99$)")
    ax.axhline(y=0.95, color="#888888", lw=0.9, linestyle=":", alpha=0.85,
               label="Acceptable stability threshold ($r_s = 0.95$)")

    # Annotate convergence regime
    if len(x) > 0 and len(y) > 0:
        sat_idx = np.where(y >= 0.99)[0]
        if len(sat_idx) > 0:
            conv_draw = x[sat_idx[0]]
            ax.scatter([conv_draw], [y[sat_idx[0]]], s=40, color="#4A7C59", zorder=5)
            ax.annotate(f"Convergence ($r_s \\geq 0.99$)\nat $n={conv_draw}$ draws",
                        xy=(conv_draw, y[sat_idx[0]]), xytext=(conv_draw + 40, y[sat_idx[0]] - 0.025),
                        arrowprops=dict(arrowstyle="->", color="#333333", lw=0.7),
                        fontsize=7, fontweight="bold", color="#333333")

    ax.set_xlabel("Dirichlet Draws ($n$)", fontsize=8, fontweight="bold")
    ax.set_ylabel("Rank Stability (Spearman $r_s$ vs Full 1,000 Draws)", fontsize=8, fontweight="bold")
    ax.set_title("Rank Stability Convergence Under Dirichlet Weight Resampling",
                 fontweight="bold", fontsize=8.5, pad=8)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.set_ylim(0.935, 1.015)
    ax.set_xlim(-10, max(x) + 30)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.subplots_adjust(left=0.15, right=0.96, top=0.90, bottom=0.14)
    save(fig, "29_convergence_analysis")


if __name__ == "__main__":
    build()
