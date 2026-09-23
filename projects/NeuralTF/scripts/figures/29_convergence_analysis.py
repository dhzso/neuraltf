"""Dirichlet draw convergence analysis.

Single-panel plot demonstrating convergence of candidate ranking stability as a function
of Dirichlet draws (k = 40) using the disjoint split-half design (two independent
half-budget draw sets, Spearman between per-gene rank vectors — no subset-vs-superset bias).
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
               label="\u03c1 = 0.99 threshold")
    ax.axhline(y=0.95, color="#888888", lw=0.9, linestyle=":", alpha=0.85,
               label="\u03c1 = 0.95 threshold")


    # Highlight the first draw count reaching the 0.99 threshold (data point only)
    if len(x) > 0 and len(y) > 0:
        sat_idx = np.where(y >= 0.99)[0]
        if len(sat_idx) > 0:
            ax.scatter([x[sat_idx[0]]], [y[sat_idx[0]]], s=40, color="#4A7C59", zorder=5)

    ax.set_xlabel("Dirichlet Draws ($n$)", fontsize=7.0)
    ax.set_ylabel("Rank Stability (Spearman $r_s$, disjoint split-half draws)", fontsize=7.0)
    # 2026-09-23: header raised ~2.5 pt + subtitle pulled 3.5 pt closer to
    # the title; legend anchor dropped 0.864 -> 0.856 so the subtitle no
    # longer cuts through the legend row (was a ~4 pt bbox overlap).
    title_block(
        fig,
        "Rank Stability Convergence Under Dirichlet Weight Resampling",
        "Split-half design: correlation of independent half-budget draw sets vs full run",
        y=0.9946, sub_y=0.9455,
    )
    h_lg, l_lg = ax.get_legend_handles_labels()
    fig.legend(handles=h_lg, labels=l_lg, fontsize=5.8, frameon=False, loc="lower center",
               bbox_to_anchor=(0.5, 0.856), ncol=3)
    ax.set_ylim(0.935, 1.015)
    ax.set_xlim(-10, max(x) + 30)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.subplots_adjust(left=0.15, right=0.96, top=0.86, bottom=0.16)
    save(fig, "29_convergence_analysis")


if __name__ == "__main__":
    build()
