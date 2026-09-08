"""Score distribution across all 9 evidence streams — genome-wide distributions.

Evaluates the dynamic range, distribution shape, and coverage of evidence scores
across all 11,675 candidates in the planarian genome.
Panel a: Genome-wide coverage (% of candidates with non-zero evidence) per stream.
Panel b: Distribution of non-zero evidence scores (violin and median) per stream.
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    df = load_all()
    streams = STREAM_COLS
    labels = [STREAM_L.get(s, s) for s in streams]
    colors = [STREAM_C.get(s, C_A) for s in streams]
    n_total = len(df)

    # Calculate coverage and extract non-zero distributions
    coverage_pct = []
    nz_scores = []
    medians = []
    for s in streams:
        if s in df.columns:
            vals = df[s].dropna().to_numpy()
            nz = vals[vals > 0]
            cov = (len(nz) / n_total) * 100
            coverage_pct.append(cov)
            nz_scores.append(nz if len(nz) > 0 else np.array([0.0]))
            medians.append(np.median(nz) if len(nz) > 0 else 0.0)
        else:
            coverage_pct.append(0.0)
            nz_scores.append(np.array([0.0]))
            medians.append(0.0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_2COL, 3.2),
                                   gridspec_kw={"width_ratios": [1, 1.4]})

    y = np.arange(len(streams))

    # Panel a: Non-zero coverage percentage
    bars = ax1.barh(y, coverage_pct, color=colors, alpha=0.85, height=0.6, edgecolor="none")
    for i, cov in enumerate(coverage_pct):
        ax1.text(cov + 1.0, i, f"{cov:.1f}%",
                 va="center", ha="left", fontsize=6.5, color="#222222")
    ax1.set_yticks(y)
    ax1.set_yticklabels(labels, fontsize=7)
    ax1.set_xlabel("Coverage (%)", fontsize=7.5)
    ax1.set_xlim(0, max(coverage_pct) * 1.18)
    ax1.set_title("Stream coverage", fontsize=8, pad=4)
    ax1.invert_yaxis()
    panel_tag(ax1, "a")

    # Panel b: Distribution of non-zero scores (violin + median dot)
    parts = ax2.violinplot(nz_scores, positions=y, vert=False, showextrema=False, widths=0.7)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(colors[i])
        pc.set_edgecolor("none")
        pc.set_alpha(0.65)

    # Add median and IQR lines
    for i, vals in enumerate(nz_scores):
        if len(vals) > 1:
            q25, med, q75 = np.percentile(vals, [25, 50, 75])
            ax2.plot([q25, q75], [i, i], color="#333333", lw=1.2)
            ax2.plot(med, i, marker="o", markersize=3.5, color="#111111", zorder=4)

    ax2.set_yticks(y)
    ax2.set_yticklabels([])  # shared labels on ax1
    ax2.set_xlabel("Score", fontsize=7.5)
    ax2.set_xlim(0, 1.05)
    ax2.set_title("Score distribution", fontsize=8, pad=4)
    ax2.invert_yaxis()
    panel_tag(ax2, "b")

    fig.tight_layout()
    save(fig, "31_score_distribution_all9")

if __name__ == "__main__":
    build()

