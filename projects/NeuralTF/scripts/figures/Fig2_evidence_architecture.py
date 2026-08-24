"""Figure 2 — Evidence architecture of prioritized candidates.

Shows WHY each candidate is highly ranked.

Panels:
  A  Annotated evidence heatmap (Top 10 × 7 streams + metadata tracks)
  B  Weighted contribution decomposition (stacked bar)
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec


STREAM_COLS = ["expression", "specificity", "reproducibility",
               "rnai", "correlation", "neural_enriched", "neural_specificity"]

WEIGHTS = {
    "expression": 0.211, "specificity": 0.105, "reproducibility": 0.158,
    "rnai": 0.158, "correlation": 0.105, "neural_enriched": 0.158,
    "neural_specificity": 0.105,
}


def _load_data():
    neural_df = load_rank_neural()
    top10 = load_top10_fixed()
    return neural_df, top10


def _order_candidates(top10, neural_df):
    """Order: Track A (by score desc), then Track B (by score desc)."""
    a = top10[top10["track"] == "A"].sort_values("composite_score", ascending=False)
    b = top10[top10["track"] == "B"].sort_values("composite_score", ascending=False)
    ordered = pd.concat([a, b])
    return ordered


def _prepare_matrix(ordered, neural_df):
    """Build a candidates × streams matrix of raw evidence values."""
    ids = ordered["gene_id"].tolist()
    matrix = []
    for gid in ids:
        row = neural_df[neural_df["gene_id"] == gid]
        if len(row) == 0:
            matrix.append([np.nan] * len(STREAM_COLS))
            continue
        vals = []
        for col in STREAM_COLS:
            v = row.iloc[0].get(col, np.nan)
            vals.append(v if pd.notna(v) else 0.0)
        matrix.append(vals)
    return np.array(matrix, dtype=float)


def _compute_contributions(matrix, weights):
    """Compute weighted contributions for each candidate × stream."""
    contributions = np.zeros_like(matrix)
    for j, stream in enumerate(STREAM_COLS):
        w = weights.get(stream, 0)
        contributions[:, j] = matrix[:, j] * w
    # Renormalize per candidate (like integrated_score does)
    for i in range(contributions.shape[0]):
        row = matrix[i]
        present = ~np.isnan(row) & (row != 0)
        if present.any():
            w_sum = sum(weights[s] for j, s in enumerate(STREAM_COLS) if present[j])
            if w_sum > 0:
                contributions[i, present] *= (1.0 / w_sum) * w_sum  # keep raw contribution
    return contributions


def fig2a_heatmap(ordered, neural_df, ax):
    """Annotated evidence heatmap: Top 10 × 7 streams with metadata tracks."""
    matrix = _prepare_matrix(ordered, neural_df)

    # Build custom colormap: white → blue gradient
    cmap = plt.cm.YlOrRd
    cmap.set_bad("#F0F0F0")

    # Plot heatmap
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1,
                   interpolation="nearest")

    # Stream labels on top
    stream_labels = [STREAM_LABELS[s] for s in STREAM_COLS]
    ax.set_xticks(np.arange(len(STREAM_COLS)))
    ax.set_xticklabels(stream_labels, rotation=45, ha="right", fontsize=5.5)
    ax.xaxis.set_label_position("top")

    # Candidate labels on left
    y_labels = []
    for _, row in ordered.iterrows():
        gid = row["gene_id"]
        name = row.get("gene_name", "")
        track = row["track"]
        prefix = "●" if track == "A" else "■"
        y_labels.append(f"{prefix} {gid}")
    ax.set_yticks(np.arange(len(ordered)))
    ax.set_yticklabels(y_labels, fontsize=5.5, fontfamily="monospace")

    # Color the y-tick labels by track
    for i, (_, row) in enumerate(ordered.iterrows()):
        color = C_TRACK_A if row["track"] == "A" else C_TRACK_B
        ax.get_yticklabels()[i].set_color(color)

    # Add cell values
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            text_color = "white" if val > 0.6 else "#333333"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=4.5, color=text_color)

    # Add vertical separator between streams
    ax.axvline(x=6.5, color="#999999", linewidth=0.5, linestyle="--")

    # Add metadata tracks on right side
    meta_x = len(STREAM_COLS) + 0.5

    # Proof status
    for i, (_, row) in enumerate(ordered.iterrows()):
        ps = str(row.get("proof_status", "")).lower()
        if "validated" in ps or "fstf" in ps.lower():
            color = PROOF_COLORS.get("prior_fstf", "#56B4E9")
        else:
            color = PROOF_COLORS.get("novel_candidate", "#E69F00")
        ax.plot(meta_x, i, "s", color=color, markersize=3, clip_on=False,
                transform=ax.get_yaxis_transform())

    # Neural specificity (if available)
    for i, (_, row) in enumerate(ordered.iterrows()):
        gid = row["gene_id"]
        neural_row = neural_df[neural_df["gene_id"] == gid]
        if len(neural_row) > 0:
            ns = neural_row.iloc[0].get("neural_specificity", np.nan)
            if pd.notna(ns) and ns > 0:
                ax.plot(meta_x + 0.5, i, "D", color=C_HIGHLIGHT, markersize=3,
                        clip_on=False, transform=ax.get_yaxis_transform())

    style_ax(ax)
    ax.set_title("Evidence heatmap", fontsize=8, fontweight="bold", pad=6, loc="left")

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.ax.set_ylabel("Evidence strength", fontsize=5)
    cbar.ax.tick_params(labelsize=5)


def fig2b_contributions(ordered, neural_df, ax):
    """Stacked bar showing weighted contribution of each stream to integrated score."""
    matrix = _prepare_matrix(ordered, neural_df)
    contributions = _compute_contributions(matrix, WEIGHTS)

    n_candidates = len(ordered)
    y = np.arange(n_candidates)[::-1]

    # Stacked horizontal bars
    left = np.zeros(n_candidates)
    for j, stream in enumerate(STREAM_COLS):
        bars = ax.barh(y, contributions[:, j], left=left, height=0.65,
                       color=STREAM_COLORS[stream], label=STREAM_LABELS[stream],
                       edgecolor="white", linewidth=0.3)
        left += contributions[:, j]

    # Candidate labels
    y_labels = []
    for _, row in ordered.iterrows():
        gid = row["gene_id"]
        track = row["track"]
        prefix = "●" if track == "A" else "■"
        y_labels.append(f"{prefix} {gid}")
    ax.set_yticks(y)
    ax.set_yticklabels(y_labels, fontsize=5.5, fontfamily="monospace")

    # Color y-tick labels by track
    for i, (_, row) in enumerate(ordered.iterrows()):
        color = C_TRACK_A if row["track"] == "A" else C_TRACK_B
        ax.get_yticklabels()[i].set_color(color)

    ax.set_xlabel("Weighted contribution to integrated score")
    style_ax(ax)
    ax.set_title("Score decomposition", fontsize=8, fontweight="bold", pad=6, loc="left")
    ax.legend(frameon=False, loc="lower right", fontsize=5, ncol=2,
              bbox_to_anchor=(1.0, -0.15))


def build():
    """Generate Figure 2 panels A–B."""
    neural_df, top10 = _load_data()
    ordered = _order_candidates(top10, neural_df)

    fig = plt.subplots(1, 1, figsize=(DOUBLE_COL, inch(70)))
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(DOUBLE_COL, inch(65)),
                                      gridspec_kw={"width_ratios": [1.2, 1.0]})

    fig2a_heatmap(ordered, neural_df, ax_a)
    fig2b_contributions(ordered, neural_df, ax_b)

    panel_letter(ax_a, "A")
    panel_letter(ax_b, "B")

    fig.tight_layout(w_pad=2.0)
    savefig(fig, "Fig2_evidence_architecture")
    return fig


if __name__ == "__main__":
    build()
