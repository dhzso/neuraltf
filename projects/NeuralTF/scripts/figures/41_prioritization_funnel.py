"""Figure 41: Prioritization Funnel — from 11,696 loci to the dual-track Top-10.

Quantitative attrition view of the NeuralTF selection cascade. Every count is
computed live from the pipeline artifacts (no hard-coded numbers):
  RUN/rank.csv        -> all annotated planarian TF loci with >=1 evidence stream
  RUN/rank_neural.csv -> neural-fate candidate universe (neural G0 mask)
  top10 CSV           -> published dual-track shortlist (5 A + 5 B)
Provenance split (RNAi-screened / previously-published FSTF / untested) is
annotated at the full-universe stage.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


def build():
    rank = load_all()
    neural = load_neural()
    top10 = load_top10()

    n_all = len(rank)
    n_neural = len(neural)
    n_top10 = len(top10)
    n_a = int((top10["track"] == "A").sum())
    n_b = int((top10["track"] == "B").sum())
    n_tested = int((rank["proof_status"] == "tested").sum())
    n_fstf = int((rank["proof_status"] == "known_fstf").sum())
    n_untested = int((rank["proof_status"] == "not_tested").sum())
    k_streams = int(rank[STREAM_COLS].notna().any(axis=1).sum())

    stages = [
        ("Annotated planarian TF loci\n(PlanMine catalog)", n_all, "#8A9AA8"),
        (f"Neural-fate candidate universe\n(neural G0 progenitor mask, 73 subclusters)",
         n_neural, C_A),
        (f"Dual-track Top-10 shortlist\n(Track A n = {n_a} screened  |  Track B n = {n_b} novel)",
         n_top10, C_B),
    ]

    fig, ax = plt.subplots(figsize=(7.08, 3.0), dpi=500)
    y = np.arange(len(stages))[::-1]
    max_n = stages[0][1]
    import math

    def width_pct(n):
        # log-scaled widths so every stage stays visible; annotated as such
        return 100.0 * math.log10(n) / math.log10(max_n)

    for yi, (label_txt, n, color) in zip(y, stages):
        ax.barh(yi, width_pct(n), left=0, color=color, height=0.58,
                edgecolor="white", lw=0.8, zorder=3)
        ax.text(width_pct(n) + 1.5, yi, f"{n:,}", va="center", ha="left",
                fontsize=9.5, fontweight="bold", color="#222222", zorder=4)
        ax.text(-1.5, yi, label_txt, va="center", ha="right", fontsize=6.4,
                color="#333333", zorder=4)

    def retention(a, b):
        return f"{b / a * 100:.2f}% retained"

    ax.annotate("", xy=(60, y[1] + 0.40), xytext=(60, y[0] - 0.40),
                arrowprops=dict(arrowstyle="->", color="#555555", lw=0.9))
    ax.text(61.5, (y[0] + y[1]) / 2, retention(n_all, n_neural) + f"  (\u2212{n_all - n_neural:,})",
            va="center", fontsize=6.2, color="#555555")
    ax.annotate("", xy=(32, y[2] + 0.40), xytext=(32, y[1] - 0.40),
                arrowprops=dict(arrowstyle="->", color="#555555", lw=0.9))
    ax.text(33.5, (y[1] + y[2]) / 2, retention(n_neural, n_top10) + f"  (\u2212{n_neural - n_top10})",
            va="center", fontsize=6.2, color="#555555")

    handles = [
        mpatches.Patch(color="#8A9AA8", label="Full annotated TF universe"),
        mpatches.Patch(color=C_A, label="Neural-fate candidates"),
        mpatches.Patch(color=C_B, label="Published dual-track shortlist"),
    ]
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, -0.03),
              frameon=False, fontsize=6.0)

    ax.text(-28, -0.42, "", fontsize=5.8)  # spacer keeps ylim breathing room

    fig.text(0.31, 0.035,
             f"Provenance at full universe: {n_tested} RNAi-screened (King 2024) · "
             f"{n_fstf} previously-published FSTF · {n_untested:,} untested",
             fontsize=5.8, color="#555555", style="italic")

    ax.set_xlim(-29, 112)
    ax.set_ylim(-0.75, len(stages) - 0.30)
    ax.set_yticks([])
    ax.set_xlabel("Share of annotated TF universe (%)  \u2014  bar width uses a log10 scale",
                  fontsize=7.0)
    ax.set_title("Prioritization Funnel: 11,696 Loci to a 10-Gene Wet-Lab Shortlist",
                 fontsize=8.5, fontweight="bold", pad=8)
    ax.spines[:].set_visible(False)
    ax.set_xticks([0, 25, 50, 75, 100])

    fig.subplots_adjust(left=0.30, right=0.97, top=0.86, bottom=0.20)
    save(fig, "41_prioritization_funnel")
    print(f"Built 41_prioritization_funnel.png ({n_all:,} -> {n_neural} -> {n_top10})")


if __name__ == "__main__":
    build()
