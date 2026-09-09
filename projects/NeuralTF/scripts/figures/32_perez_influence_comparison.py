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

    data = [neural_cls, other_cls, unclass]
    med_n = np.median(neural_cls)
    med_o = np.median(other_cls)
    med_u = np.median(unclass)

    labels = [
        f"Neural Lineage\n(n={len(neural_cls):,}, med={med_n:.2f})",
        f"Other Lineages\n(n={len(other_cls):,}, med={med_o:.2f})",
        f"Unclassified\n(n={len(unclass):,}, med={med_u:.2f})",
    ]
    colors = [C_A, C_B, "#B0BEC5"]

    # Boxplot
    bp = ax.boxplot(data, patch_artist=True, widths=0.48, zorder=4,
                    medianprops=dict(color="#111111", lw=1.4),
                    boxprops=dict(lw=0.8),
                    whiskerprops=dict(lw=0.8, color="#444444"),
                    capprops=dict(lw=0.8, color="#444444"),
                    flierprops=dict(marker="none"))

    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col)
        patch.set_alpha(0.70)
        patch.set_edgecolor("#333333")

    # Overlay jitter for all 3 groups
    np.random.seed(42)
    j1 = np.random.normal(0, 0.045, size=len(neural_cls))
    ax.scatter(1 + j1, neural_cls, s=12, color=C_A, alpha=0.35, edgecolors="none", zorder=3)

    j2 = np.random.normal(0, 0.045, size=len(other_cls))
    ax.scatter(2 + j2, other_cls, s=12, color=C_B, alpha=0.35, edgecolors="none", zorder=3)

    # Subsample n=300 points for unclassified group to show consistent point representation
    unclass_sub = np.random.choice(unclass, size=300, replace=False)
    j3 = np.random.normal(0, 0.045, size=len(unclass_sub))
    ax.scatter(3 + j3, unclass_sub, s=8, color="#788896", alpha=0.20, edgecolors="none", zorder=3)

    # Statistical tests
    stat1, pval1 = mannwhitneyu(neural_cls, other_cls, alternative="two-sided")
    stat2, pval2 = mannwhitneyu(neural_cls, unclass, alternative="two-sided")
    p1_str = r"P < 10^{-15}" if pval1 < 1e-15 else f"P = {pval1:.1e}"
    p2_str = r"P < 10^{-15}" if pval2 < 1e-15 else f"P = {pval2:.1e}"

    # Significance bracket 1: Neural vs Other (1 to 2)
    y_bar1 = 1.04
    h = 0.025
    ax.plot([1, 1, 2, 2], [y_bar1, y_bar1 + h, y_bar1 + h, y_bar1], color="#333333", lw=0.7)
    ax.text(1.5, y_bar1 + h + 0.015, f"Mann–Whitney: ${p1_str}$",
            ha="center", va="bottom", fontsize=6.2, color="#222222")

    # Significance bracket 2: Neural vs Unclassified (1 to 3)
    y_bar2 = 1.16
    ax.plot([1, 1, 3, 3], [y_bar2, y_bar2 + h, y_bar2 + h, y_bar2], color="#333333", lw=0.7)
    ax.text(2.0, y_bar2 + h + 0.015, f"Mann–Whitney: ${p2_str}$",
            ha="center", va="bottom", fontsize=6.2, color="#222222")

    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(labels, fontsize=6.8)
    ax.set_ylabel("Integrated Evidence Score", fontsize=7.0, fontweight="bold")
    ax.set_title("Evidence Score Stratification Across Single-Cell Lineage Classes",
                 fontsize=8.0, pad=8, fontweight="bold")
    ax.set_ylim(-0.02, 1.30)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.subplots_adjust(left=0.14, right=0.96, top=0.88, bottom=0.16)
    save(fig, "32_perez_influence_comparison")


if __name__ == "__main__":
    build()
