"""Figure 3 — Ranking robustness under weight uncertainty.

Shows how candidate rankings behave under three weighting methods.

Panels:
  A  Bump chart: Top 10 across fixed, centered Dirichlet, uniform Dirichlet
  B  Rank uncertainty from weight sensitivity draws
  C  Candidate robustness classification
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def _load_data():
    fixed = load_top10_fixed()
    centered = load_top10_centered()
    uniform = load_top10_uniform()
    draws = load_weight_sensitivity_draws()
    top10_ch = load_weight_sensitivity_top10()
    return fixed, centered, uniform, draws, top10_ch


def _build_rank_table(fixed, centered, uniform):
    """Build a table of candidates × methods with their ranks."""
    # Collect all unique candidates from the three top-10 lists
    all_ids = set()
    for df in [fixed, centered, uniform]:
        if "gene_id" in df.columns:
            all_ids.update(df["gene_id"].tolist())

    # For each method, assign ranks by composite_score
    records = []
    for gid in all_ids:
        rec = {"gene_id": gid}

        # Fixed: rank among fixed candidates
        if "gene_id" in fixed.columns:
            f_row = fixed[fixed["gene_id"] == gid]
            if len(f_row) > 0:
                rec["fixed_rank"] = f_row.iloc[0].get("rank", np.nan)
                rec["fixed_score"] = f_row.iloc[0].get("composite_score",
                                          f_row.iloc[0].get("integrated_score", np.nan))
                rec["track_fixed"] = f_row.iloc[0].get("track", "")

        # Centered
        if "gene_id" in centered.columns:
            c_row = centered[centered["gene_id"] == gid]
            if len(c_row) > 0:
                rec["centered_rank"] = c_row.iloc[0].get("rank", np.nan)
                rec["centered_score"] = c_row.iloc[0].get("composite_score",
                                            c_row.iloc[0].get("dirichlet_median_score", np.nan))

        # Uniform
        if "gene_id" in uniform.columns:
            u_row = uniform[uniform["gene_id"] == gid]
            if len(u_row) > 0:
                rec["uniform_rank"] = u_row.iloc[0].get("rank", np.nan)
                rec["uniform_score"] = u_row.iloc[0].get("composite_score",
                                           u_row.iloc[0].get("uniform_median_score", np.nan))

        records.append(rec)

    return pd.DataFrame(records)


def fig3a_bump_chart(rank_table, fixed, ax):
    """Parallel-rank (bump) chart showing rank shifts across three methods."""
    methods = ["Fixed", "Centered", "uniform"]
    rank_cols = ["fixed_rank", "centered_rank", "uniform_rank"]

    # Only show candidates that appear in at least 2 methods
    present = rank_table.dropna(subset=rank_cols, thresh=2)
    if present.empty:
        ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes,
                ha="center", va="center", fontsize=8)
        return

    # Track assignment from fixed
    track_map = {}
    for _, row in fixed.iterrows():
        track_map[row["gene_id"]] = row.get("track", "")

    x = np.arange(len(methods))
    for _, row in present.iterrows():
        gid = row["gene_id"]
        ranks = [row.get(c, np.nan) for c in rank_cols]
        track = track_map.get(gid, "")

        # Skip candidates with no valid ranks
        valid = [(i, r) for i, r in enumerate(ranks) if pd.notna(r)]
        if len(valid) < 2:
            continue

        x_valid = [v[0] for v in valid]
        y_valid = [v[1] for v in valid]

        color = C_TRACK_A if track == "A" else C_TRACK_B
        alpha = 1.0 if gid in fixed["gene_id"].values else 0.5

        ax.plot(x_valid, y_valid, "-o", color=color, alpha=alpha,
                markersize=4, linewidth=1.0, zorder=3)

        # Label at the rightmost point
        if x_valid:
            last_x = x_valid[-1]
            last_y = y_valid[-1]
            label = gene_label(fixed, gid) if "gene_name" in fixed.columns else gid
            ax.annotate(label, (last_x, last_y), fontsize=4.5,
                       xytext=(5, 0), textcoords="offset points",
                       va="center", ha="left", color=color, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(["Fixed\nweight", "Dirichlet\ncentered", "Dirichlet\nuniform"],
                       fontsize=6)
    ax.set_ylabel("Rank in Top 10")
    ax.invert_yaxis()
    ax.set_ylim(0.5, 10.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    style_ax(ax, title="")
    ax.set_title("Rank stability across methods", fontsize=8, fontweight="bold",
                 pad=6, loc="left")

    # Legend
    legend_elements = [
        Line2D([0], [0], color=C_TRACK_A, marker="o", label="Track A", linewidth=1),
        Line2D([0], [0], color=C_TRACK_B, marker="o", label="Track B", linewidth=1),
    ]
    ax.legend(handles=legend_elements, frameon=False, fontsize=5, loc="lower right")


def fig3b_rank_uncertainty(draws, fixed, ax):
    """Show rank uncertainty intervals from weight sensitivity draws."""
    # Get top 10 from fixed
    top_ids = fixed["gene_id"].tolist()[:10]

    # For each candidate, compute rank distribution across draws
    gene_draws = draws[draws["gene_id"].isin(top_ids)]

    if gene_draws.empty:
        ax.text(0.5, 0.5, "No sensitivity data", transform=ax.transAxes,
                ha="center", va="center", fontsize=8)
        return

    # Compute statistics per candidate
    stats = gene_draws.groupby("gene_id")["rank"].agg(
        median_rank="median",
        q25=lambda x: np.percentile(x, 2.5),
        q75=lambda x: np.percentile(x, 75),
        min_rank="min",
        max_rank="max",
        in_top10_frac=lambda x: (x <= 10).mean()
    ).reset_index()

    # Order by median rank
    stats = stats.sort_values("median_rank")
    y = np.arange(len(stats))[::-1]

    # Track colors
    track_map = dict(zip(fixed["gene_id"], fixed.get("track", [""] * len(fixed))))
    colors = [C_TRACK_A if track_map.get(gid, "") == "A" else C_TRACK_B
              for gid in stats["gene_id"]]

    # Uncertainty intervals (25th-75th percentile)
    ax.barh(y, stats["q75"] - stats["q25"], left=stats["q25"],
            height=0.5, color=colors, alpha=0.3, edgecolor="none")

    # Full range whiskers
    for i, (_, row) in enumerate(stats.iterrows()):
        ax.plot([row["min_rank"], row["max_rank"]], [y[i], y[i]],
                color=colors[i], linewidth=0.5, alpha=0.5)

    # Median rank
    ax.scatter(stats["median_rank"], y, color=colors, s=12, zorder=5,
               edgecolors="white", linewidth=0.3)

    # Labels
    y_labels = [gene_label(fixed, gid) for gid in stats["gene_id"]]
    ax.set_yticks(y)
    ax.set_yticklabels(y_labels, fontsize=5.5)

    # Annotate P(Top 10)
    for i, (_, row) in enumerate(stats.iterrows()):
        frac = row["in_top10_frac"]
        ax.text(ax.get_xlim()[1] + 0.3, y[i], f"{frac:.0%}",
                fontsize=5, va="center", ha="left",
                color="#0072B2" if frac >= 0.8 else "#999999")

    ax.set_xlabel("Rank position")
    ax.set_title("Rank uncertainty (1000 weight draws)", fontsize=8,
                 fontweight="bold", pad=6, loc="left")
    style_ax(ax)

    # Add P(Top 10) label
    ax.text(1.02, 1.02, "P(Top10)", transform=ax.transAxes, fontsize=5,
            fontweight="bold", va="bottom", ha="left")


def fig3c_robustness_classification(top10_ch, fixed, ax):
    """Classify candidates by robustness: stable, moderate, sensitive."""
    if top10_ch.empty:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=8)
        return

    df = top10_ch.copy()
    # Classify
    def classify(row):
        frac = row.get("frac_draws_in_top10", 0)
        if frac >= 0.8:
            return "Stable"
        elif frac >= 0.5:
            return "Moderate"
        else:
            return "Sensitive"

    df["class"] = df.apply(classify, axis=1)

    class_colors = {"Stable": "#0072B2", "Moderate": "#E69F00", "Sensitive": "#D55E00"}
    class_order = ["Stable", "Moderate", "Sensitive"]

    for cls in class_order:
        subset = df[df["class"] == cls]
        if subset.empty:
            continue
        y = np.arange(len(subset))
        bars = ax.barh(y + class_order.index(cls) * 0.1, subset["frac_draws_in_top10"],
                       height=0.6, color=class_colors[cls], alpha=0.7, label=cls)

        # Labels
        for i, (_, row) in enumerate(subset.iterrows()):
            label = row.get("gene_name", row["gene_id"])
            ax.text(0.02, y[i] + class_order.index(cls) * 0.1,
                    f"  {label}", fontsize=5, va="center", ha="left")

    ax.set_xlabel("Fraction of draws in Top 10")
    ax.set_yticks([])
    ax.set_xlim(0, 1.05)
    style_ax(ax, title="")
    ax.set_title("Robustness classification", fontsize=8, fontweight="bold",
                 pad=6, loc="left")
    ax.legend(frameon=False, fontsize=5, loc="lower right")


def build():
    """Generate Figure 3 panels A–C."""
    fixed, centered, uniform, draws, top10_ch = _load_data()
    rank_table = _build_rank_table(fixed, centered, uniform)

    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE_COL, inch(60)),
                             gridspec_kw={"width_ratios": [1.2, 1.2, 1.0]})

    fig3a_bump_chart(rank_table, fixed, axes[0])
    fig3b_rank_uncertainty(draws, fixed, axes[1])
    fig3c_robustness_classification(top10_ch, fixed, axes[2])

    for i, ax in enumerate(axes):
        panel_letter(ax, chr(65 + i))

    fig.tight_layout(w_pad=1.5)
    savefig(fig, "Fig3_ranking_robustness")
    return fig


if __name__ == "__main__":
    build()
