"""Comprehensive evidence stream heatmap for all 134 neural-filtered transcription factors.

Single-panel, publication-quality 500 DPI figure integrating:
- Biological track categorization (Track A: RNAi-validated, Track B: Novel candidates, Prior FSTFs)
- 9 continuous and specialized evidence streams ordered to prevent visual discreteness
- Aligned candidate-specific integrated evidence scores
- Clean standardized gene symbols for all 134 candidates
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


def build():
    neural = load_neural()

    # 1. Partition into tracks and sort by integrated_score descending
    t_a = neural[neural["proof_status"] == "known_rnai_validated"].sort_values("integrated_score", ascending=False)
    t_b = neural[neural["proof_status"] == "novel_candidate"].sort_values("integrated_score", ascending=False)
    t_f = neural[neural["proof_status"] == "prior_fstf_not_tested"].sort_values("integrated_score", ascending=False)
    df = pd.concat([t_a, t_b, t_f], ignore_index=True)

    # 2. Arrange streams: continuous dense streams first, then specialized and functional
    streams = [
        "expression",
        "specificity",
        "neural_specificity",
        "reproducibility",
        "correlation",
        "perez_lineage",
        "perez_influence",
        "rnai",
        "neural_enriched",
    ]
    stream_labels = [
        "Expression",
        "Specificity",
        "Neural Spec.",
        "Reproducibility",
        "Correlation",
        "Lineage Assoc.",
        "Network Centrality",
        "RNAi Validation",
        "Neural Filter",
    ]

    mat = df[streams].fillna(0).values
    n_rows = len(df)
    n_cols = len(streams)

    # Soft, publication-grade colormap: neutral pearl gray -> pale steel -> ocean blue -> midnight navy
    cmap_colors = ["#F6F8FA", "#D5E3EE", "#97BDD7", "#5792BF", "#2C6492", "#16324F"]
    custom_cmap = mcolors.LinearSegmentedColormap.from_list("pub_blues", cmap_colors)

    # Figure dimensions: 7.2 x 15.5 inches for 134 rows at 500 DPI
    fig = plt.figure(figsize=(7.2, 15.5), dpi=500)
    gs = fig.add_gridspec(1, 3, width_ratios=[0.02, 0.45, 0.47], wspace=0.18)
    ax_track = fig.add_subplot(gs[0, 0])
    ax_main = fig.add_subplot(gs[0, 1])
    ax_score = fig.add_subplot(gs[0, 2], sharey=ax_main)

    C_FSTF = "#7C786E"
    colors = [C_A] * len(t_a) + [C_B] * len(t_b) + [C_FSTF] * len(t_f)
    div1 = len(t_a) - 0.5
    div2 = len(t_a) + len(t_b) - 0.5

    # 1. Track strip
    for i, c in enumerate(colors):
        ax_track.add_patch(plt.Rectangle((0, i - 0.5), 1, 1, color=c, ec="none"))
    ax_track.set_xlim(0, 1)
    ax_track.set_ylim(n_rows - 0.5, -0.5)
    ax_track.axhline(div1, color="white", lw=1.5)
    ax_track.axhline(div2, color="white", lw=1.5)
    ax_track.axis("off")

    # 2. Main heatmap
    im = ax_main.imshow(mat, aspect="auto", cmap=custom_cmap, vmin=0, vmax=1, interpolation="nearest")
    ax_main.set_xticks(range(n_cols))
    ax_main.set_xticklabels(stream_labels, rotation=42, ha="left", fontsize=6.8, fontweight="bold")
    ax_main.xaxis.tick_top()
    ax_main.xaxis.set_label_position("top")
    ax_main.set_ylim(n_rows - 0.5, -0.5)
    ax_main.axhline(div1, color="white", lw=1.5)
    ax_main.axhline(div2, color="white", lw=1.5)
    ax_main.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    for s in ax_main.spines.values():
        s.set_visible(False)

    # 3. Score bars with clean gene labels
    scores = df["integrated_score"].values
    y_pos = np.arange(n_rows)
    gene_labels = [clean_gene_symbol(r.gene_name, r.gene_id) for _, r in df.iterrows()]

    ax_score.barh(y_pos, scores, height=0.74, color=colors, edgecolor="none", alpha=0.92)
    ax_score.set_xlim(0, 1.15)
    ax_score.set_ylim(n_rows - 0.5, -0.5)
    ax_score.axhline(div1, color="#CCCCCC", lw=0.8, ls="--")
    ax_score.axhline(div2, color="#CCCCCC", lw=0.8, ls="--")
    ax_score.set_xlabel("Integrated Evidence Score", fontsize=7.2, fontweight="bold")
    ax_score.xaxis.tick_top()
    ax_score.xaxis.set_label_position("top")
    ax_score.tick_params(axis="x", labelsize=6.5)

    ax_score.set_yticks(y_pos)
    ax_score.set_yticklabels(gene_labels, fontsize=5.2)
    ax_score.tick_params(axis="y", pad=6, length=0)

    ax_score.spines["top"].set_visible(True)
    ax_score.spines["top"].set_linewidth(0.6)
    ax_score.spines["right"].set_visible(False)
    ax_score.spines["bottom"].set_visible(False)
    ax_score.spines["left"].set_visible(True)
    ax_score.spines["left"].set_color("#D0D0D0")
    ax_score.spines["left"].set_linewidth(0.5)

    # Annotate numeric score at end of each bar
    for y, s in enumerate(scores):
        ax_score.text(s + 0.02, y, f"{s:.2f}", va="center", ha="left", fontsize=4.6, color="#333333")

    # Layout adjustment
    fig.subplots_adjust(left=0.06, right=0.95, top=0.93, bottom=0.045)

    # Horizontal Colorbar below heatmap
    cbar_ax = fig.add_axes([0.10, 0.02, 0.35, 0.009])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Evidence Stream Score", fontsize=6.5, fontweight="bold")
    cbar.ax.tick_params(labelsize=6)
    cbar.outline.set_linewidth(0.5)

    # Track Legend below score bars
    leg_handles = [
        Patch(facecolor=C_A, label=f"Track A: RNAi-validated (n={len(t_a)})"),
        Patch(facecolor=C_B, label=f"Track B: Novel candidates (n={len(t_b)})"),
        Patch(facecolor=C_FSTF, label=f"Prior FSTF: Untested (n={len(t_f)})"),
    ]
    fig.legend(
        handles=leg_handles,
        loc="lower right",
        bbox_to_anchor=(0.95, 0.015),
        ncol=1,
        frameon=False,
        fontsize=6.8,
    )

    # Group label text on the far left margin
    fig.text(0.018, 0.65, f"Track A: Validated (n={len(t_a)})", rotation=90, va="center", ha="center", fontsize=7.5, fontweight="bold", color=C_A)
    fig.text(0.018, 0.28, f"Track B: Novel (n={len(t_b)})", rotation=90, va="center", ha="center", fontsize=7.5, fontweight="bold", color=C_B)
    fig.text(0.018, 0.08, f"Prior FSTF (n={len(t_f)})", rotation=90, va="center", ha="center", fontsize=7.5, fontweight="bold", color=C_FSTF)

    fig.suptitle(
        "Evidence Stream Integration Across All 134 Prioritized Neural Transcription Factors",
        fontweight="bold",
        fontsize=8.8,
        y=0.985,
    )

    save(fig, "04_evidence_heatmap_neural")


if __name__ == "__main__":
    build()
