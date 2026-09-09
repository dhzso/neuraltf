"""Comprehensive evidence stream heatmap for all 134 neural-filtered transcription factors.

Publication-quality single-panel 500 DPI figure integrating:
- Biological track categorization (Track A: RNAi-validated, Track B: Novel candidates, Prior FSTFs)
- 7 informative continuous & specialized evidence streams ordered to maintain visual coherence
- Distinguishes missing/unmeasured values with neutral gray fill instead of false zeros
- Integrated evidence score representation with clean scale and no text clutter
- Standardized, unique gene symbols for all 134 candidates (with isoform disambiguation)
"""
from __future__ import annotations
import sys, re
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


def get_unique_gene_label(row, all_df):
    """Return clean gene symbol, adding isoform suffix (.1, .2) only for multi-isoform loci."""
    base = clean_gene_symbol(row.get("gene_name"), row.get("gene_id"))
    gid = str(row.get("gene_id", ""))
    
    # Check if this base symbol appears multiple times across all candidates
    same_base = all_df[all_df["base_symbol"] == base]
    if len(same_base) > 1:
        # Extract isoform suffix from gene_id (e.g. _0_1 -> .1, _0_2 -> .2)
        m = re.search(r"_0_([0-9]+)$", gid)
        if m:
            return f"{base}.{m.group(1)}"
    return base


def build():
    neural = load_neural()

    # Add temporary base symbol for isoform disambiguation
    neural["base_symbol"] = [clean_gene_symbol(r.gene_name, r.gene_id) for _, r in neural.iterrows()]

    # 1. Partition into biological tracks, sorted by integrated_score descending within track
    t_a = neural[neural["proof_status"] == "known_rnai_validated"].sort_values("integrated_score", ascending=False)
    t_b = neural[neural["proof_status"] == "novel_candidate"].sort_values("integrated_score", ascending=False)
    t_f = neural[neural["proof_status"] == "prior_fstf_not_tested"].sort_values("integrated_score", ascending=False)
    df = pd.concat([t_a, t_b, t_f], ignore_index=True)

    # Clean unique gene labels
    gene_labels = [get_unique_gene_label(r, neural) for _, r in df.iterrows()]

    # 2. Select 7 active evidence streams:
    # Continuous high-density streams first, then specialized network/correlation streams
    streams = [
        "expression",
        "specificity",
        "neural_specificity",
        "reproducibility",
        "perez_lineage",
        "perez_influence",
        "correlation",
    ]
    stream_labels = [
        "Bulk Expression",
        "Bulk Specificity",
        "Neural Spec. (scRNA)",
        "Reproducibility",
        "Lineage Association",
        "Network Centrality",
        "Marker Co-expression",
    ]

    # Extract raw data and mask NaNs (unmeasured/not applicable)
    raw_mat = df[streams].values.astype(float)
    masked_mat = np.ma.masked_invalid(raw_mat)
    n_rows = len(df)
    n_cols = len(streams)

    # Colormap: soft elegant gradient from ice white-blue to deep ocean navy
    cmap_colors = ["#EDF4F9", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C", "#08306B"]
    base_cmap = mcolors.LinearSegmentedColormap.from_list("pub_blues", cmap_colors)
    # Neutral light gray for missing/unmeasured values
    base_cmap.set_bad(color="#ECEFF1")

    # Dimensions: 8.0 x 12.0 inches at 500 DPI (crisp, balanced publication format)
    fig = plt.figure(figsize=(8.0, 12.0), dpi=500)
    gs = fig.add_gridspec(
        1, 4,
        width_ratios=[0.022, 0.50, 0.20, 0.26],
        wspace=0.10,
        left=0.06,
        right=0.96,
        top=0.87,
        bottom=0.05
    )
    ax_track = fig.add_subplot(gs[0, 0])
    ax_main = fig.add_subplot(gs[0, 1])
    ax_score = fig.add_subplot(gs[0, 2], sharey=ax_main)
    ax_labels = fig.add_subplot(gs[0, 3], sharey=ax_main)

    C_FSTF = "#706E65"
    colors = [C_A] * len(t_a) + [C_B] * len(t_b) + [C_FSTF] * len(t_f)
    div1 = len(t_a) - 0.5
    div2 = len(t_a) + len(t_b) - 0.5

    # 1. Track strip (far left)
    for i, c in enumerate(colors):
        ax_track.add_patch(plt.Rectangle((0, i - 0.5), 1, 1, color=c, ec="none"))
    ax_track.set_xlim(0, 1)
    ax_track.set_ylim(n_rows - 0.5, -0.5)
    ax_track.axhline(div1, color="white", lw=2.0)
    ax_track.axhline(div2, color="white", lw=2.0)
    ax_track.axis("off")

    # 2. Main heatmap
    im = ax_main.imshow(masked_mat, aspect="auto", cmap=base_cmap, vmin=0, vmax=1, interpolation="nearest")
    ax_main.set_xticks(range(n_cols))
    ax_main.set_xticklabels(stream_labels, rotation=35, ha="left", fontsize=6.8, fontweight="bold")
    ax_main.xaxis.tick_top()
    ax_main.xaxis.set_label_position("top")
    ax_main.tick_params(axis="x", pad=3)
    ax_main.set_ylim(n_rows - 0.5, -0.5)
    ax_main.axhline(div1, color="white", lw=2.0)
    ax_main.axhline(div2, color="white", lw=2.0)
    # Subtle separator between dense streams (0..4) and specialized streams (5..6)
    ax_main.axvline(4.5, color="#B0BEC5", lw=0.8, ls=":")
    ax_main.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    for s in ax_main.spines.values():
        s.set_visible(False)

    # 3. Integrated Score bars (clean, elegant, no row-by-row text numbers)
    scores = df["integrated_score"].values
    y_pos = np.arange(n_rows)

    ax_score.barh(y_pos, scores, height=0.72, color=colors, edgecolor="none", alpha=0.88)
    ax_score.set_xlim(0, 1.05)
    ax_score.set_ylim(n_rows - 0.5, -0.5)
    ax_score.axhline(div1, color="#D0D0D0", lw=1.0, ls="--")
    ax_score.axhline(div2, color="#D0D0D0", lw=1.0, ls="--")
    # Reference threshold lines at 0.70 and 0.85
    ax_score.axvline(0.70, color="#9E9E9E", lw=0.6, ls=":")
    ax_score.set_xlabel("Integrated Score", fontsize=7.2, fontweight="bold", labelpad=6)
    ax_score.xaxis.tick_top()
    ax_score.xaxis.set_label_position("top")
    ax_score.set_xticks([0.0, 0.5, 1.0])
    ax_score.tick_params(axis="x", labelsize=6.5, pad=3)
    ax_score.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    ax_score.spines["top"].set_visible(True)
    ax_score.spines["top"].set_linewidth(0.6)
    ax_score.spines["right"].set_visible(False)
    ax_score.spines["bottom"].set_visible(False)
    ax_score.spines["left"].set_visible(False)

    # 4. Clean Gene Labels on the right
    ax_labels.set_ylim(n_rows - 0.5, -0.5)
    ax_labels.set_xlim(0, 1)
    ax_labels.axis("off")
    for y, (name, col) in enumerate(zip(gene_labels, colors)):
        # Key candidates bolded slightly
        weight = "bold" if y < 5 or (len(t_a) <= y < len(t_a) + 5) else "normal"
        ax_labels.text(0.04, y, name, va="center", ha="left", fontsize=4.8, fontweight=weight, color="#222222")

    # Colorbar and Missing Value Legend below heatmap
    cbar_ax = fig.add_axes([0.14, 0.022, 0.28, 0.009])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Normalized Evidence Stream Score", fontsize=6.5, fontweight="bold")
    cbar.set_ticks([0.0, 0.5, 1.0])
    cbar.ax.tick_params(labelsize=6)
    cbar.outline.set_linewidth(0.5)

    # Missing / Not evaluated swatch
    fig.patches.append(plt.Rectangle((0.44, 0.022), 0.015, 0.009, transform=fig.transFigure,
                                     facecolor="#ECEFF1", edgecolor="#CCCCCC", lw=0.5, clip_on=False))
    fig.text(0.46, 0.025, "Not evaluated / N/A", fontsize=6.0, va="center", color="#555555")

    # Track Legend below score bars
    leg_handles = [
        Patch(facecolor=C_A, label=f"Track A: RNAi-validated (n={len(t_a)})"),
        Patch(facecolor=C_B, label=f"Track B: Novel candidates (n={len(t_b)})"),
        Patch(facecolor=C_FSTF, label=f"Prior FSTF: Untested (n={len(t_f)})"),
    ]
    fig.legend(
        handles=leg_handles,
        loc="lower right",
        bbox_to_anchor=(0.96, 0.015),
        ncol=1,
        frameon=False,
        fontsize=6.5,
    )

    # Track labels on far left margin
    fig.text(0.02, 0.62, f"Track A: Validated (n={len(t_a)})", rotation=90, va="center", ha="center", fontsize=7.5, fontweight="bold", color=C_A)
    fig.text(0.02, 0.25, f"Track B: Novel (n={len(t_b)})", rotation=90, va="center", ha="center", fontsize=7.5, fontweight="bold", color=C_B)
    fig.text(0.02, 0.075, f"Prior FSTF (n={len(t_f)})", rotation=90, va="center", ha="center", fontsize=7.5, fontweight="bold", color=C_FSTF)

    fig.suptitle(
        "Multi-Stream Evidence Landscape Across All 134 Prioritized Neural Transcription Factors",
        fontweight="bold",
        fontsize=8.8,
        y=0.985,
    )

    save(fig, "04_evidence_heatmap_neural")


if __name__ == "__main__":
    build()

