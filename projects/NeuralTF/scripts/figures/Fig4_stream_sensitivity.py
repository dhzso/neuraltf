"""Figure 4 — Sensitivity of evidence streams.

Shows the effect of removing each evidence stream on candidate prioritization.

Panels:
  A  Global impact: rank correlation change when each stream is removed
  B  Candidate-specific sensitivity heatmap: rank change per candidate × stream
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


STREAM_COLS = ["expression", "specificity", "reproducibility",
               "rnai", "correlation", "neural_enriched", "neural_specificity"]

WEIGHTS = {
    "expression": 0.211, "specificity": 0.105, "reproducibility": 0.158,
    "rnai": 0.158, "correlation": 0.105, "neural_enriched": 0.158,
    "neural_specificity": 0.105,
}

DEFAULT_WEIGHTS = np.array([WEIGHTS[s] for s in STREAM_COLS])


def _load_data():
    neural_df = load_rank_neural()
    top10 = load_top10_fixed()
    return neural_df, top10


def _compute_scores_with_stream_removed(df, stream_idx):
    """Compute integrated score with one stream removed."""
    scores = np.zeros(len(df))
    for i, (_, row) in enumerate(df.iterrows()):
        vals = np.array([row.get(s, 0) if pd.notna(row.get(s, np.nan)) else np.nan
                        for s in STREAM_COLS])
        weights = DEFAULT_WEIGHTS.copy()
        weights[stream_idx] = 0  # Remove this stream

        present = ~np.isnan(vals) & (vals != 0)
        if present.any():
            w = weights[present]
            v = vals[present]
            scores[i] = np.sum(w * v) / np.sum(w)
        else:
            scores[i] = 0
    return scores


def _compute_ablation_results(neural_df, top10):
    """For each stream removal, compute rank changes for all 99 candidates."""
    # Baseline scores (with all streams)
    baseline_scores = np.zeros(len(neural_df))
    for i, (_, row) in enumerate(neural_df.iterrows()):
        vals = np.array([row.get(s, 0) if pd.notna(row.get(s, np.nan)) else np.nan
                        for s in STREAM_COLS])
        present = ~np.isnan(vals) & (vals != 0)
        if present.any():
            w = DEFAULT_WEIGHTS[present]
            v = vals[present]
            baseline_scores[i] = np.sum(w * v) / np.sum(w)

    baseline_ranks = pd.Series(baseline_scores).rank(ascending=False).values

    results = []
    for j, stream in enumerate(STREAM_COLS):
        ablated_scores = _compute_scores_with_stream_removed(neural_df, j)
        ablated_ranks = pd.Series(ablated_scores).rank(ascending=False).values

        rank_changes = ablated_ranks - baseline_ranks  # positive = rank got worse

        for i, (_, row) in enumerate(neural_df.iterrows()):
            results.append({
                "gene_id": row["gene_id"],
                "stream_removed": stream,
                "rank_change": rank_changes[i],
                "score_change": ablated_scores[i] - baseline_scores[i],
                "baseline_rank": baseline_ranks[i],
                "ablated_rank": ablated_ranks[i],
                "in_top10_baseline": baseline_ranks[i] <= 10,
                "in_top10_ablated": ablated_ranks[i] <= 10,
                "displaced_from_top10": (baseline_ranks[i] <= 10) and (ablated_ranks[i] > 10),
            })

    return pd.DataFrame(results)


def fig4a_global_impact(ablation_df, ax):
    """Global impact: Spearman correlation and median |rank change| per stream."""
    streams = STREAM_COLS
    n_streams = len(streams)

    # Aggregate per stream
    stats = []
    for stream in streams:
        sub = ablation_df[ablation_df["stream_removed"] == stream]
        median_abs_change = sub["rank_change"].abs().median()
        mean_abs_change = sub["rank_change"].abs().mean()
        n_displaced = sub["displaced_from_top10"].sum()
        n_top10_changes = ((sub["in_top10_baseline"] != sub["in_top10_ablated"])).sum()
        stats.append({
            "stream": stream,
            "median_abs_rank_change": median_abs_change,
            "mean_abs_rank_change": mean_abs_change,
            "n_displaced": n_displaced,
            "n_top10_changes": n_top10_changes,
        })

    stats_df = pd.DataFrame(stats)

    # Lollipop plot: median absolute rank change
    y = np.arange(n_streams)[::-1]
    x = stats_df["median_abs_rank_change"].values
    colors = [STREAM_COLORS[s] for s in streams]

    ax.barh(y, x, height=0.5, color=colors, alpha=0.7, edgecolor="white", linewidth=0.3)
    ax.scatter(x, y, color=colors, s=20, zorder=5, edgecolors="white", linewidth=0.5)

    ax.set_yticks(y)
    ax.set_yticklabels([STREAM_LABELS[s] for s in streams], fontsize=6)
    ax.set_xlabel("Median |rank change| (all 99 candidates)")

    # Annotate n_top10_changes
    for i, (xi, n) in enumerate(zip(x, stats_df["n_top10_changes"])):
        ax.text(xi + 0.1, y[i], f"ΔTop10={n}", fontsize=5, va="center", ha="left",
                color="#D55E00" if n > 0 else "#999999")

    style_ax(ax, title="")
    ax.set_title("Global stream impact", fontsize=8, fontweight="bold", pad=6, loc="left")


def fig4b_candidate_sensitivity(ablation_df, top10, ax):
    """Heatmap: rank change per candidate × stream removed.
    Only show candidates that appear in any top-10 baseline."""
    # Filter to top-10 candidates
    top10_ids = set(top10["gene_id"].tolist())
    sub = ablation_df[ablation_df["gene_id"].isin(top10_ids)]

    if sub.empty:
        ax.text(0.5, 0.5, "No top-10 candidates in ablation data",
                transform=ax.transAxes, ha="center", va="center", fontsize=8)
        return

    # Pivot: candidates × streams
    pivot = sub.pivot_table(index="gene_id", columns="stream_removed",
                            values="rank_change", aggfunc="first")

    # Reorder streams
    pivot = pivot.reindex(columns=[s for s in STREAM_COLS if s in pivot.columns])

    # Order candidates by track
    track_map = dict(zip(top10["gene_id"], top10.get("track", [""] * len(top10))))
    pivot["_track"] = pivot.index.map(lambda x: track_map.get(x, ""))
    pivot = pivot.sort_values(["_track", "gene_id"], ascending=[True, True])
    pivot = pivot.drop(columns=["_track"])

    # Diverging colormap centered on 0
    vmax = max(abs(pivot.values.min()), abs(pivot.values.max()), 1)
    cmap = plt.cm.RdBu_r
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    im = ax.imshow(pivot.values, aspect="auto", cmap=cmap, norm=norm,
                   interpolation="nearest")

    # Labels
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels([STREAM_LABELS[s] for s in pivot.columns],
                       rotation=45, ha="right", fontsize=5.5)
    ax.xaxis.set_label_position("top")

    y_labels = [gene_label(top10, gid) for gid in pivot.index]
    ax.set_yticks(np.arange(len(pivot)))
    ax.set_yticklabels(y_labels, fontsize=5.5)

    # Color y-tick labels by track
    for i, gid in enumerate(pivot.index):
        color = C_TRACK_A if track_map.get(gid, "") == "A" else C_TRACK_B
        ax.get_yticklabels()[i].set_color(color)

    # Add cell values
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.values[i, j]
            if not np.isnan(val):
                text_color = "white" if abs(val) > vmax * 0.6 else "#333333"
                ax.text(j, i, f"{val:.0f}", ha="center", va="center",
                        fontsize=4.5, color=text_color)

    style_ax(ax, title="")
    ax.set_title("Candidate sensitivity", fontsize=8, fontweight="bold", pad=6, loc="left")

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.ax.set_ylabel("Rank change (− = improves)", fontsize=5)
    cbar.ax.tick_params(labelsize=5)


def build():
    """Generate Figure 4 panels A–B."""
    neural_df, top10 = _load_data()
    ablation_df = _compute_ablation_results(neural_df, top10)

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL, inch(65)),
                             gridspec_kw={"width_ratios": [1.0, 1.2]})

    fig4a_global_impact(ablation_df, axes[0])
    fig4b_candidate_sensitivity(ablation_df, top10, axes[1])

    panel_letter(axes[0], "A")
    panel_letter(axes[1], "B")

    fig.tight_layout(w_pad=1.5)
    savefig(fig, "Fig4_stream_sensitivity")
    return fig


if __name__ == "__main__":
    build()
