"""Lineage-resolved evidence enrichment — integrated score by single-cell lineage class.

Single-panel boxplot with jitter evaluating whether single-cell lineage classification
from Perez et al. enriches for high integrated evidence scores.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def build():
    rank = load_all()

    if "perez_lineage" not in rank.columns or not rank["perez_lineage"].notna().any():
        raise FileNotFoundError("rank.csv carries no perez_lineage values.")

    fig, ax = plt.subplots(figsize=(W_15COL, 3.6), dpi=500)

    neural_cls = rank[rank["perez_lineage"] == 1.0]["integrated_score"].dropna().values
    other_cls = rank[rank["perez_lineage"] == 0.5]["integrated_score"].dropna().values
    unclass = rank[rank["perez_lineage"] == 0.0]["integrated_score"].dropna().values

    groups = [
        ("Neural Lineage\n(Perez et al.)", neural_cls, C_A),
        ("Other Lineages\n(Non-neural)", other_cls, C_B),
        ("Unclassified\n(Transcriptome)", unclass, "#D0D7DE"),
    ]

    data = [g[1] for g in groups]
    labels = [f"{g[0]}\n(n={len(g[1]):,})" for g in groups]
    colors = [g[2] for g in groups]

    # Boxplot
    bp = ax.boxplot(data, patch_artist=True, widths=0.52,
                    medianprops=dict(color="#111111", lw=1.3),
                    boxprops=dict(lw=0.8),
                    whiskerprops=dict(lw=0.8, color="#555555"),
                    capprops=dict(lw=0.8, color="#555555"),
                    flierprops=dict(marker=".", markersize=2.5, alpha=0.25, color="#888888"))

    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col)
        patch.set_alpha(0.75)
        patch.set_edgecolor("#333333")

    # Overlay subtle jitter for small groups
    np.random.seed(42)
    for i, vals in enumerate([neural_cls, other_cls]):
        jitter = np.random.normal(0, 0.04, size=len(vals))
        ax.scatter(i + 1 + jitter, vals, s=12, color=colors[i], alpha=0.45,
                   edgecolors="none", zorder=3)

    # Statistical test between Neural and Other lineages
    stat, pval = mannwhitneyu(neural_cls, other_cls, alternative="two-sided")
    p_str = "P < 10^{-15}" if pval < 1e-15 else f"P = {pval:.1e}"

    # Significance bracket between group 1 and group 2
    y_bar = max(np.percentile(neural_cls, 95), np.percentile(other_cls, 95)) + 0.12
    h = 0.02
    ax.plot([1, 1, 2, 2], [y_bar, y_bar + h, y_bar + h, y_bar], color="#222222", lw=0.8)
    ax.text(1.5, y_bar + h + 0.015, f"Mann–Whitney U: {p_str}",
            ha="center", va="bottom", fontsize=6.2, color="#222222")

    ax.set_xticklabels(labels, fontsize=6.8)
    ax.set_ylabel("Integrated Evidence Score", fontsize=7.0)
    ax.set_title("Evidence Score Stratification Across Single-Cell Lineage Classes",
                 fontsize=8.0, pad=6)
    ax.set_ylim(-0.02, 1.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.subplots_adjust(left=0.14, right=0.96, top=0.90, bottom=0.16)
    save(fig, "32_perez_influence_comparison")


if __name__ == "__main__":
    build()
