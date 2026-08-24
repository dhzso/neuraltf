"""Figure 1 — Candidate prioritization landscape.

Shows how 249 TF candidates are progressively reduced to the final
prioritized Top 10.

Panels:
  A  Filtering flow (alluvial-style horizontal bars)
  B  Score ECDF with highlighted populations
  C  Ranked landscape with final candidates marked
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D


def _load_data():
    all_df = load_rank_all()
    neural_df = load_rank_neural()
    top10 = load_top10_fixed()
    uniform_rank = load_uniform_full_rank()

    # Merge uniform scores into all_df for comparison
    if "uniform_median_score" in uniform_rank.columns:
        all_df = all_df.merge(
            uniform_rank[["gene_id", "uniform_median_score"]],
            on="gene_id", how="left"
        )
    return all_df, neural_df, top10


def _count_populations(all_df, neural_df, top10):
    n_all = len(all_df)
    # Expression filter: p <= 0.05
    expr_pass = (all_df["expression"] > 0).sum() if "expression" in all_df.columns else n_all
    n_neural = len(neural_df)
    track_a = top10[top10["track"] == "A"]
    track_b = top10[top10["track"] == "B"]
    n_a = len(track_a)
    n_b = len(track_b)
    return {
        "all": n_all,
        "expression": expr_pass,
        "neural": n_neural,
        "track_a": n_a,
        "track_b": n_b,
        "top10": n_a + n_b,
    }


def fig1a_filtering_flow(counts, ax):
    """Horizontal bar chart showing population at each filtering stage."""
    stages = [
        "All TF candidates",
        "Expression-supported",
        "Neural-filtered",
        "Track A (RNAi-validated)",
        "Track B (novel candidates)",
        "Final Top 10",
    ]
    values = [
        counts["all"],
        counts["expression"],
        counts["neural"],
        counts["track_a"],
        counts["track_b"],
        counts["top10"],
    ]
    colors = [C_ALL249, C_NEURAL, C_CENTERED, C_TRACK_A, C_TRACK_B, C_HIGHLIGHT]

    y = np.arange(len(stages))[::-1]
    bars = ax.barh(y, values, color=colors, edgecolor="white", linewidth=0.5, height=0.6)

    for bar, val in zip(bars, values):
        pct = val / counts["all"] * 100
        ax.text(bar.get_width() + counts["all"] * 0.01, bar.get_y() + bar.get_height() / 2,
                f"{val} ({pct:.0f}%)", va="center", fontsize=6, color="#333333")

    ax.set_yticks(y)
    ax.set_yticklabels(stages, fontsize=6)
    ax.set_xlim(0, counts["all"] * 1.25)
    ax.set_xlabel("Number of candidates")
    style_ax(ax, title="")
    ax.set_title("Filtering stages", fontsize=8, fontweight="bold", pad=6, loc="left")

    # Add connecting lines for flow effect
    for i in range(len(stages) - 1):
        y_top = y[i] + 0.3
        y_bot = y[i + 1] + 0.3
        w_top = values[i]
        w_bot = values[i + 1]
        ax.plot([w_top, w_top], [y_top, y_top + 0.05], color="#999999", linewidth=0.3)
        ax.plot([w_bot, w_bot], [y_bot - 0.05, y_bot], color="#999999", linewidth=0.3)
        ax.plot([w_top, w_bot], [y_top + 0.05, y_bot - 0.05],
                color="#999999", linewidth=0.3, linestyle="--", alpha=0.5)


def fig1b_score_ecdf(all_df, neural_df, top10, ax):
    """ECDF of integrated scores for all 249, highlighting neural and top-10."""
    score_col = "integrated_score"
    if score_col not in all_df.columns:
        # Try alternative
        for c in ["fixed_weight_score", "score"]:
            if c in all_df.columns:
                score_col = c
                break

    # Sort for ECDF
    all_scores = np.sort(all_df[score_col].dropna().values)
    ecdf_y = np.arange(1, len(all_scores) + 1) / len(all_scores)

    ax.step(all_scores, ecdf_y, where="post", color=C_ALL249, linewidth=1.0, label="All 249 TFs")

    # Neural subset
    if score_col in neural_df.columns:
        neural_scores = np.sort(neural_df[score_col].dropna().values)
        ecdf_n = np.arange(1, len(neural_scores) + 1) / len(neural_scores)
        ax.step(neural_scores, ecdf_n, where="post", color=C_CENTERED, linewidth=1.0,
                label="99 neural-filtered")

    # Top 10 markers
    top_ids = top10["gene_id"].tolist()
    for track, color, label in [("A", C_TRACK_A, "Track A"), ("B", C_TRACK_B, "Track B")]:
        track_ids = top10[top10["track"] == track]["gene_id"].tolist()
        track_df = neural_df[neural_df["gene_id"].isin(track_ids)]
        if score_col in track_df.columns and len(track_df) > 0:
            vals = track_df[score_col].values
            for v in vals:
                # Find ECDF value from all_scores
                rank = np.searchsorted(all_scores, v, side="right") / len(all_scores)
                ax.plot(v, rank, "o", color=color, markersize=4, zorder=5)

    ax.set_xlabel("Integrated evidence score")
    ax.set_ylabel("Cumulative fraction")
    style_ax(ax, title="")
    ax.set_title("Score distributions", fontsize=8, fontweight="bold", pad=6, loc="left")
    ax.legend(frameon=False, loc="lower right", fontsize=5)

    # Add population annotations
    y_ann = 0.15
    for label_text, color, count in [
        ("All 249", C_ALL249, len(all_df)),
        ("Neural 99", C_CENTERED, len(neural_df)),
        ("Track A", C_TRACK_A, len(top10[top10["track"] == "A"])),
        ("Track B", C_TRACK_B, len(top10[top10["track"] == "B"])),
    ]:
        ax.text(0.97, y_ann, f"{label_text}: n={count}", transform=ax.transAxes,
                fontsize=5, ha="right", va="top", color=color, fontweight="bold")
        y_ann -= 0.06


def fig1c_ranked_landscape(all_df, neural_df, top10, ax):
    """Scatter of rank vs integrated score for all candidates,
    with final Top 10 highlighted."""
    score_col = "integrated_score"
    if score_col not in all_df.columns:
        for c in ["fixed_weight_score", "score"]:
            if c in all_df.columns:
                score_col = c
                break

    # Rank all candidates by score
    all_sorted = all_df.sort_values(score_col, ascending=False).reset_index(drop=True)
    all_sorted["rank_pos"] = np.arange(1, len(all_sorted) + 1)

    # Background: all 249
    ax.scatter(all_sorted["rank_pos"], all_sorted[score_col],
               s=4, c=C_ALL249, alpha=0.4, edgecolors="none", label="All 249 TFs")

    # Neural-filtered
    neural_ids = set(neural_df["gene_id"].tolist())
    neural_mask = all_sorted["gene_id"].isin(neural_ids)
    ax.scatter(all_sorted.loc[neural_mask, "rank_pos"],
               all_sorted.loc[neural_mask, score_col],
               s=6, c=C_CENTERED, alpha=0.5, edgecolors="none", label="Neural 99")

    # Top 10
    for track, color, marker in [("A", C_TRACK_A, "o"), ("B", C_TRACK_B, "s")]:
        track_ids = set(top10[top10["track"] == track]["gene_id"].tolist())
        track_mask = all_sorted["gene_id"].isin(track_ids)
        if track_mask.any():
            ax.scatter(all_sorted.loc[track_mask, "rank_pos"],
                       all_sorted.loc[track_mask, score_col],
                       s=30, c=color, marker=marker, edgecolors="white",
                       linewidth=0.5, zorder=5, label=f"Track {'A' if track == 'A' else 'B'}")

    # Annotate Top 10
    top10_ids = set(top10["gene_id"].tolist())
    annotate_df = all_sorted[all_sorted["gene_id"].isin(top10_ids)]
    for _, row in annotate_df.iterrows():
        label = gene_label(all_df, row["gene_id"])
        ax.annotate(label, (row["rank_pos"], row[score_col]),
                    fontsize=4.5, ha="center", va="bottom",
                    xytext=(0, 4), textcoords="offset points",
                    fontweight="bold", color="#333333")

    ax.set_xlabel("Rank (by integrated score)")
    ax.set_ylabel("Integrated evidence score")
    style_ax(ax, title="")
    ax.set_title("Ranked landscape", fontsize=8, fontweight="bold", pad=6, loc="left")
    ax.legend(frameon=False, loc="lower left", fontsize=5, markerscale=1.5)


def build():
    """Generate Figure 1 panels A–C."""
    all_df, neural_df, top10 = _load_data()
    counts = _count_populations(all_df, neural_df, top10)

    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE_COL, inch(60)),
                             gridspec_kw={"width_ratios": [1.0, 1.0, 1.3]})

    fig1a_filtering_flow(counts, axes[0])
    fig1b_score_ecdf(all_df, neural_df, top10, axes[1])
    fig1c_ranked_landscape(all_df, neural_df, top10, axes[2])

    for i, ax in enumerate(axes):
        panel_letter(ax, chr(65 + i))

    fig.tight_layout(w_pad=1.5)
    savefig(fig, "Fig1_candidate_landscape")
    return fig


if __name__ == "__main__":
    build()
