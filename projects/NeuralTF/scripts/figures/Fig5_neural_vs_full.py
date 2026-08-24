"""Figure 5 — Neural filtering vs full candidate universe.

Demonstrates that neural-filtered candidates are robustly recovered
in the full 249-candidate universe.

Panels:
  A  Rank-rank comparison (99-neural vs 249-wide)
  B  Score distribution comparison with overlap
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
from matplotlib.patches import Patch


def _load_data():
    all_df = load_rank_all()
    neural_df = load_rank_neural()
    top10 = load_top10_fixed()
    uniform_99 = load_uniform_full_rank()
    uniform_249 = load_uniform_all249_full_rank()
    return all_df, neural_df, top10, uniform_99, uniform_249


def fig5a_rank_comparison(all_df, neural_df, top10, uniform_99, uniform_249, ax):
    """Rank-rank scatter: 99-neural uniform rank vs 249-wide uniform rank."""
    if uniform_99.empty or uniform_249.empty:
        ax.text(0.5, 0.5, "Missing uniform rank data", transform=ax.transAxes,
                ha="center", va="center", fontsize=8)
        return

    # Both use uniform_median_score
    score_col_99 = "uniform_median_score"
    score_col_249 = "uniform_median_score"

    if score_col_99 not in uniform_99.columns or score_col_249 not in uniform_249.columns:
        ax.text(0.5, 0.5, "Missing score columns", transform=ax.transAxes,
                ha="center", va="center", fontsize=8)
        return

    # Rank within each scope
    u99 = uniform_99.copy()
    u99["rank_99"] = u99[score_col_99].rank(ascending=False)

    u249 = uniform_249.copy()
    u249["rank_249"] = u249[score_col_249].rank(ascending=False)

    # Merge on gene_id
    merged = u99[["gene_id", "rank_99"]].merge(
        u249[["gene_id", "rank_249"]], on="gene_id", how="inner"
    )

    if merged.empty:
        ax.text(0.5, 0.5, "No overlapping candidates", transform=ax.transAxes,
                ha="center", va="center", fontsize=8)
        return

    # Top 10 identification
    top10_ids = set(top10["gene_id"].tolist())
    track_map = dict(zip(top10["gene_id"], top10.get("track", [""] * len(top10))))

    # Background
    non_top10 = merged[~merged["gene_id"].isin(top10_ids)]
    ax.scatter(non_top10["rank_99"], non_top10["rank_249"],
               s=6, c=C_NEURAL, alpha=0.3, edgecolors="none", label="Other neural candidates")

    # Top 10
    for track, color in [("A", C_TRACK_A), ("B", C_TRACK_B)]:
        t_ids = [gid for gid, t in track_map.items() if t == track]
        sub = merged[merged["gene_id"].isin(t_ids)]
        if not sub.empty:
            ax.scatter(sub["rank_99"], sub["rank_249"],
                       s=30, c=color, edgecolors="white", linewidth=0.5,
                       zorder=5, label=f"Track {'A' if track == 'A' else 'B'}")

    # Annotate top 10
    for _, row in merged[merged["gene_id"].isin(top10_ids)].iterrows():
        label = gene_label(all_df, row["gene_id"])
        ax.annotate(label, (row["rank_99"], row["rank_249"]),
                    fontsize=4.5, ha="center", va="bottom",
                    xytext=(0, 4), textcoords="offset points",
                    fontweight="bold", color="#333333")

    # Identity line
    lims = [0, max(merged["rank_99"].max(), merged["rank_249"].max()) * 1.05]
    ax.plot(lims, lims, "--", color="#999999", linewidth=0.5, zorder=1)

    # Spearman correlation
    from scipy.stats import spearmanr
    rho, pval = spearmanr(merged["rank_99"], merged["rank_249"])
    ax.text(0.95, 0.05, f"Spearman ρ = {rho:.3f}\np = {pval:.1e}",
            transform=ax.transAxes, fontsize=6, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    ax.set_xlabel("Rank (99 neural, uniform Dirichlet)")
    ax.set_ylabel("Rank (249 wide, uniform Dirichlet)")
    ax.legend(frameon=False, fontsize=5, loc="upper left")
    style_ax(ax, title="")
    ax.set_title("Rank recovery across scopes", fontsize=8, fontweight="bold",
                 pad=6, loc="left")


def fig5b_score_comparison(all_df, neural_df, top10, uniform_99, uniform_249, ax):
    """Score distributions: 99 neural vs 249 wide, with top-10 highlighted."""
    score_col = "integrated_score"
    if score_col not in all_df.columns:
        for c in ["fixed_weight_score", "score"]:
            if c in all_df.columns:
                score_col = c
                break

    # Neural scores
    neural_scores = neural_df[score_col].dropna().values if score_col in neural_df.columns else []
    all_scores = all_df[score_col].dropna().values if score_col in all_df.columns else []

    # Top 10 scores
    top10_ids = set(top10["gene_id"].tolist())
    track_map = dict(zip(top10["gene_id"], top10.get("track", [""] * len(top10))))

    # ECDF for all 249
    if len(all_scores) > 0:
        sorted_all = np.sort(all_scores)
        ecdf_all = np.arange(1, len(sorted_all) + 1) / len(sorted_all)
        ax.step(sorted_all, ecdf_all, where="post", color=C_ALL249, linewidth=1.0,
                label=f"All 249 (n={len(all_scores)})")

    # ECDF for 99 neural
    if len(neural_scores) > 0:
        sorted_neural = np.sort(neural_scores)
        ecdf_neural = np.arange(1, len(sorted_neural) + 1) / len(sorted_neural)
        ax.step(sorted_neural, ecdf_neural, where="post", color=C_CENTERED, linewidth=1.0,
                label=f"Neural 99 (n={len(neural_scores)})")

    # Top 10 markers
    for track, color, marker in [("A", C_TRACK_A, "o"), ("B", C_TRACK_B, "s")]:
        track_ids = [gid for gid, t in track_map.items() if t == track]
        track_df = neural_df[neural_df["gene_id"].isin(track_ids)]
        if score_col in track_df.columns and not track_df.empty:
            for v in track_df[score_col].values:
                if len(all_scores) > 0:
                    rank = np.searchsorted(sorted_all, v, side="right") / len(sorted_all)
                    ax.plot(v, rank, marker=marker, color=color, markersize=5,
                            zorder=5, markeredgecolor="white", markeredgewidth=0.3)

    ax.set_xlabel("Integrated evidence score")
    ax.set_ylabel("Cumulative fraction")
    style_ax(ax, title="")
    ax.set_title("Score distributions by scope", fontsize=8, fontweight="bold",
                 pad=6, loc="left")
    ax.legend(frameon=False, fontsize=5, loc="lower right")

    # Add KS test
    if len(all_scores) > 0 and len(neural_scores) > 0:
        from scipy.stats import ks_2samp
        ks_stat, ks_p = ks_2samp(all_scores, neural_scores)
        ax.text(0.95, 0.35, f"KS = {ks_stat:.3f}\np = {ks_p:.2e}",
                transform=ax.transAxes, fontsize=6, ha="right", va="bottom",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))


def build():
    """Generate Figure 5 panels A–B."""
    all_df, neural_df, top10, uniform_99, uniform_249 = _load_data()

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL, inch(55)))

    fig5a_rank_comparison(all_df, neural_df, top10, uniform_99, uniform_249, axes[0])
    fig5b_score_comparison(all_df, neural_df, top10, uniform_99, uniform_249, axes[1])

    panel_letter(axes[0], "A")
    panel_letter(axes[1], "B")

    fig.tight_layout(w_pad=1.5)
    savefig(fig, "Fig5_neural_filtering_scope")
    return fig


if __name__ == "__main__":
    build()
