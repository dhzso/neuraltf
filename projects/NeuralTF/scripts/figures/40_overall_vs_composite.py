"""Figure 40: Integrated Score vs Composite Score (bonus impact).

Scatter plot showing how the domain/ortholog/GO bonuses shift candidate
rankings. Top-10 published shortlist highlighted by track.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def build():
    all_df = pd.read_csv(RUN / "rank.csv")
    fixed = pd.read_csv(RES / "fixed_full_rank.csv")
    top10 = load_top10()

    # Merge composite_score into the rank frame
    merged = all_df.merge(
        fixed[["gene_id", "composite_score", "bonus_total", "rank"]],
        on="gene_id", how="left", suffixes=("", "_fixed")
    )

    # Rank by integrated_score
    merged["int_rank"] = merged["integrated_score"].rank(ascending=False, method="first").astype(int)

    top10_ids = set(top10["gene_id"].values)
    is_top10 = merged["gene_id"].isin(top10_ids)
    top10_df = merged[is_top10].merge(top10[["gene_id", "track"]], on="gene_id", how="left")
    bg_df = merged[~is_top10]

    fig, ax = plt.subplots(figsize=(4.4, 4.0), dpi=500)

    # Background: all candidates
    ax.scatter(bg_df["integrated_score"], bg_df["composite_score"],
               s=4, c="#D0D0D0", alpha=0.4, edgecolors="none", zorder=1,
               label=f"Other candidates (n = {len(bg_df):,})")

    # Diagonal reference
    ax.plot([0, 1], [0, 1], color="#999999", lw=0.5, ls="--", zorder=0,
            label="Composite = integrated (y = x)")

    # Top-10 by track
    track_colors = {"A": C_A, "B": C_HL}
    for track, grp in top10_df.groupby("track"):
        c = track_colors.get(track, "#888888")
        ax.scatter(grp["integrated_score"], grp["composite_score"],
                   s=45, c=c, edgecolors="white", linewidths=0.8, zorder=5,
                   label=f"Track {track} (n = {len(grp)})")

    # Annotate top-10 in a clean right-side column with leader lines (no overlap)
    ann = top10_df.sort_values("composite_score", ascending=False).reset_index(drop=True)
    y_top, y_bot = 1.05, 0.60
    ys = np.linspace(y_top, y_bot, len(ann))
    for yi, (_, row) in zip(ys, ann.iterrows()):
        name = label(top10_df, row["gene_id"])
        c = track_colors.get(row.get("track", "?"), "#333")
        ax.annotate(name,
                    xy=(row["integrated_score"], row["composite_score"]),
                    xytext=(0.50, yi), textcoords="data",
                    fontsize=5.2, fontweight="bold", color=c,
                    ha="left", va="center", zorder=6,
                    arrowprops=dict(arrowstyle="-", color=c, lw=0.4, alpha=0.7))

    ax.set_xlabel("Integrated evidence score", fontsize=7, fontweight="bold")
    ax.set_ylabel("Composite score (+ bonuses)", fontsize=7, fontweight="bold")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.08)
    ax.legend(fontsize=5.5, loc="upper left", frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Annotation: mean bonus shift
    delta = (top10_df["composite_score"] - top10_df["integrated_score"]).mean()
    title_block(
        fig,
        "Impact of Additive Bonuses on Candidate Scores",
        f"Mean bonus $\\Delta$ = +{delta:.3f} (composite \u2212 integrated); dashed = identity $y = x$; $N$ = {len(merged):,} loci plotted",
    )

    fig.subplots_adjust(left=0.15, right=0.95, top=0.86, bottom=0.14)
    save(fig, "40_overall_vs_composite")
    print("Built 40_overall_vs_composite.png")


if __name__ == "__main__":
    build()
