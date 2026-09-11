"""Score distributions across all 9 evidence streams for the planarian genome.

Single-panel horizontal violin plot showing empirical density, median, and interquartile range (IQR)
of non-zero evidence scores for each evidence stream across all 11,675 candidates.
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
    df = load_all()
    streams = STREAM_COLS
    labels = [STREAM_L.get(s, s) for s in streams]
    colors = [STREAM_C.get(s, C_A) for s in streams]

    nz_scores = []
    medians = []
    counts = []

    for s in streams:
        if s in df.columns:
            vals = df[s].dropna().to_numpy()
            nz = vals[vals > 0]
            nz_scores.append(nz if len(nz) > 0 else np.array([0.0]))
            medians.append(np.median(nz) if len(nz) > 0 else 0.0)
            counts.append(len(nz))
        else:
            nz_scores.append(np.array([0.0]))
            medians.append(0.0)
            counts.append(0)

    fig, ax = plt.subplots(figsize=(W_15COL, 4.2), dpi=500)
    y = np.arange(len(streams))

    # Horizontal violin distributions
    parts = ax.violinplot(nz_scores, positions=y, orientation="horizontal", showextrema=False, widths=0.72)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(colors[i])
        pc.set_edgecolor("none")
        pc.set_alpha(0.65)

    # Median and IQR bars
    for i, vals in enumerate(nz_scores):
        if len(vals) > 1:
            q25, med, q75 = np.percentile(vals, [25, 50, 75])
            ax.plot([q25, q75], [i, i], color="#222222", lw=1.5, zorder=3)
            ax.plot(med, i, marker="o", markersize=4.0, color="#111111", zorder=4)
        elif len(vals) == 1:
            ax.plot(vals[0], i, marker="o", markersize=4.0, color="#111111", zorder=4)

        # Annotation of median score and non-zero count
        ax.text(1.04, i, f"Med={medians[i]:.2f} (n={counts[i]:,})",
                va="center", ha="left", fontsize=6.0, color="#333333")

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.set_xlabel("Normalized Evidence Score", fontsize=7.0)
    ax.set_xlim(-0.02, 1.35)
    ax.set_title("Empirical Score Distributions Across 11 Integrated Evidence Streams",
                 fontsize=8.0, pad=6)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.subplots_adjust(left=0.22, right=0.96, top=0.90, bottom=0.12)
    save(fig, "31_score_distribution_all9")


if __name__ == "__main__":
    build()
