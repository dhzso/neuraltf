"""Figure 38: Leave-One-Atlas-Out Robustness.

Bar chart of per-stream shortlist Jaccard overlap and Spearman correlation
when each atlas/evidence stream is excluded in turn.
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
    summary = pd.read_csv(RES / "loo_atlas_summary.csv")

    # Sort by shortlist_overlap (ascending) so most disruptive is at top
    summary = summary.sort_values("shortlist_overlap", ascending=True)

    streams = summary["excluded_stream"].values
    labels = [STREAM_L.get(s, s) for s in streams]
    jaccard = summary["shortlist_jaccard"].values
    spearman = summary["spearman_correlation"].values
    n_new = summary["shortlist_new_in"].apply(lambda x: len(eval(x)) if isinstance(x, str) else 0).values
    n_drop = summary["shortlist_dropped"].apply(lambda x: len(eval(x)) if isinstance(x, str) else 0).values

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.08, 3.5), dpi=500,
                                     gridspec_kw={"width_ratios": [1.2, 1]})

    # Panel A: Shortlist Jaccard overlap
    colors = [C_A if j >= 0.8 else C_HL if j < 0.6 else "#8B7355" for j in jaccard]
    bars = ax1.barh(range(len(labels)), jaccard, color=colors, edgecolor="white", lw=0.5, height=0.7)
    ax1.set_yticks(range(len(labels)))
    ax1.set_yticklabels(labels, fontsize=6)
    ax1.set_xlabel("Shortlist Jaccard overlap  (in / out = gained / lost)", fontsize=7, fontweight="bold")
    ax1.set_ylabel("Excluded evidence stream", fontsize=7, fontweight="bold")
    ax1.set_xlim(0, 1.42)
    ax1.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax1.axvline(x=1.0, color="#999999", lw=0.5, ls=":", zorder=0)
    for i, (j, nn, nd) in enumerate(zip(jaccard, n_new, n_drop)):
        ax1.text(j + 0.03, i, f"{j:.2f}  ({nn} in / {nd} out)", va="center", fontsize=5.4)
    ax1.set_title("A  Shortlist stability", fontsize=8, fontweight="bold", loc="left")
    ax1.invert_yaxis()
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: Spearman correlation of full rank
    colors2 = [C_A if s >= 0.95 else "#8B7355" if s >= 0.85 else C_HL for s in spearman]
    ax2.barh(range(len(labels)), spearman, color=colors2, edgecolor="white", lw=0.5, height=0.7)
    ax2.set_yticks(range(len(labels)))
    ax2.set_yticklabels(labels, fontsize=6)
    ax2.set_xlabel("Spearman ρ (full rank)", fontsize=7, fontweight="bold")
    ax2.set_xlim(0.70, 1.15)
    ax2.set_xticks([0.70, 0.80, 0.90, 1.00])
    ax2.axvline(x=1.0, color="#999999", lw=0.5, ls=":", zorder=0)
    for i, s in enumerate(spearman):
        ax2.text(s + 0.008, i, f"{s:.3f}", va="center", fontsize=5.4)
    ax2.set_title("B  Rank correlation", fontsize=8, fontweight="bold", loc="left")
    ax2.invert_yaxis()
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    sub_bot = title_block(
        fig,
        "Leave-One-Stream-Out Robustness",
        "Each stream excluded once; shortlist and full-rank stability recomputed",
    )
    # 2026-09-23: single-row legends placed in a dedicated band between the
    # subtitle and the panel titles (the old 3-row legends anchored at axes-frac
    # 1.015 hung ~0.75 in above each axes and collided with title/subtitle).
    _pt = 1.0 / (72.0 * fig.get_size_inches()[1])
    _ax_b, _ax_t = 0.14, 0.74
    _leg_bottom = sub_bot - (6.0 + 10.25) * _pt  # 6 pt gap + 1-row legend (fs 5.0)
    _leg_y = (_leg_bottom - _ax_b) / (_ax_t - _ax_b)  # axes-fraction anchor
    from matplotlib.patches import Patch
    ax1.legend(handles=[Patch(facecolor=C_A, label="≥ 0.80"),
                        Patch(facecolor="#8B7355", label="0.60–0.80"),
                        Patch(facecolor=C_HL, label="< 0.60")],
               frameon=False, fontsize=5.0, loc="lower right", ncol=3,
               bbox_to_anchor=(1.0, _leg_y))
    ax2.legend(handles=[Patch(facecolor=C_A, label="≥ 0.95"),
                        Patch(facecolor="#8B7355", label="0.85–0.95"),
                        Patch(facecolor=C_HL, label="< 0.85")],
               frameon=False, fontsize=5.0, loc="lower right", ncol=3,
               bbox_to_anchor=(1.0, _leg_y))
    fig.subplots_adjust(left=0.17, right=0.98, top=_ax_t, bottom=_ax_b, wspace=0.62)
    save(fig, "38_loo_robustness")
    print("Built 38_loo_robustness.png")


if __name__ == "__main__":
    build()
