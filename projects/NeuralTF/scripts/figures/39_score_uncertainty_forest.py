"""Figure 39: Score-Uncertainty Forest Plot.

Horizontal error bars showing the 95% weight-perturbation band for each
top-10 candidate under Dirichlet-centered weight uncertainty. Candidates
are ordered by integrated score (highest at top).
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
    ci = pd.read_csv(RES / "bootstrap_scores_ci.csv")
    top10 = load_top10()

    # Merge to get track info and limit to top-10
    top10_ids = set(top10["gene_id"].values)
    ci_top = ci[ci["gene_id"].isin(top10_ids)].copy()
    ci_top = ci_top.merge(top10[["gene_id", "track"]], on="gene_id", how="left")
    ci_top = ci_top.sort_values("integrated_score", ascending=True)  # bottom-to-top

    # Build labels
    labels = []
    for _, row in ci_top.iterrows():
        name = label(ci_top, row["gene_id"])
        track = row.get("track", "?")
        labels.append(f"{name}  (Track {track})")

    mean = ci_top["centered_mean"].values
    lo = ci_top["centered_ci_95_lo"].values
    hi = ci_top["centered_ci_95_hi"].values
    point = ci_top["integrated_score"].values
    band_width = hi - lo

    fig, ax = plt.subplots(figsize=(3.5, 3.8), dpi=500)

    # Color by track
    track_colors = {"A": C_A, "B": C_HL}
    y_positions = np.arange(len(labels))

    for i, (y, m, l, h, pt, tr) in enumerate(
            zip(y_positions, mean, lo, hi, point, ci_top["track"].values)):
        c = track_colors.get(tr, "#888888")
        # Band
        ax.plot([l, h], [y, y], color=c, lw=1.8, solid_capstyle="round")
        # Point estimate
        ax.scatter(pt, y, color=c, s=28, zorder=5, edgecolors="white", linewidths=0.5)
        # Band width annotation
        ax.text(h + 0.008, y, f"[{l:.3f}, {h:.3f}]  Δ={h-l:.3f}",
                va="center", fontsize=5.2, color="#555555")

    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlabel("Integrated Evidence Score", fontsize=7, fontweight="bold")
    ax.set_xlim(0.62, 1.17)
    ax.set_title("Weight-Perturbation Bands (Dirichlet k=40)", fontsize=8,
                 fontweight="bold", pad=16)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend above the axes (outside the data area — no collision)
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=C_A, lw=2, label="Track A (RNAi-screened)"),
        Line2D([0], [0], color=C_HL, lw=2, label="Track B (Unscreened)"),
    ]
    ax.legend(handles=legend_elements, fontsize=5.5, loc="lower center",
              bbox_to_anchor=(0.5, 1.01), ncol=2, frameon=False,
              handletextpad=0.4, columnspacing=1.4)

    fig.subplots_adjust(left=0.32, right=0.92, top=0.93, bottom=0.10)
    save(fig, "39_score_uncertainty_forest")
    print("Built 39_score_uncertainty_forest.png")


if __name__ == "__main__":
    build()
