"""Candidate rank stability trajectories across Dirichlet prior weighting schemes.

Single-panel bump chart tracking rank conservation from fixed weights to centered Dirichlet (k=40)
and uniform Dirichlet (alpha=1) for top prioritized candidates.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


def build():
    f10 = load_top10()
    cf = load_centered_full()
    uf = load_uniform_full()
    neural = load_neural()

    fig, ax = plt.subplots(figsize=(W_15COL, 4.6), dpi=500)

    methods = ["Fixed Weight", "Dirichlet Centered\n(k=40)", "Dirichlet Uniform\n(\u03b1=1)"]
    x_pos = [0, 1, 2]

    left_labels = []
    right_labels = []

    for _, r in f10.iterrows():
        gid = r["gene_id"]
        nm = clean_gene_symbol(r.get("gene_name", ""), gid)
        tr = r["track"]
        color = C_A if tr == "A" else C_B
        ls = "-" if tr == "A" else "--"

        r_fix = neural[neural["gene_id"] == gid].index[0] + 1
        c_sub = cf[cf["gene_id"].isin(neural["gene_id"])].reset_index(drop=True)
        r_cen = c_sub[c_sub["gene_id"] == gid].index[0] + 1
        u_sub = uf[uf["gene_id"].isin(neural["gene_id"])].reset_index(drop=True)
        r_uni = u_sub[u_sub["gene_id"] == gid].index[0] + 1

        ranks = [r_fix, r_cen, r_uni]
        ax.plot(x_pos, ranks, marker="o", color=color, linestyle=ls, lw=1.5,
                markersize=5, alpha=0.9, zorder=4)
        left_labels.append((r_fix, f"{nm} (#{r_fix})", color))
        right_labels.append((r_uni, f"#{r_uni} {nm}", color))

    # Declutter left labels
    left_labels.sort(key=lambda x: x[0])
    adj_left = []
    last_y = -999
    for y, txt, c in left_labels:
        cur_y = max(y, last_y + 2.2)
        adj_left.append((cur_y, txt, c))
        last_y = cur_y

    for y, txt, c in adj_left:
        ax.text(-0.06, y, txt, ha="right", va="center", fontsize=6.2, color=c)

    # Declutter right labels
    right_labels.sort(key=lambda x: x[0])
    adj_right = []
    last_y = -999
    for y, txt, c in right_labels:
        cur_y = max(y, last_y + 2.2)
        adj_right.append((cur_y, txt, c))
        last_y = cur_y

    for y, txt, c in adj_right:
        ax.text(2.06, y, txt, ha="left", va="center", fontsize=6.2, color=c)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(methods, fontsize=6.8)
    ax.set_ylabel("Neural Candidate Rank (1–134)", fontsize=7.0)
    ax.set_ylim(-3, 72)
    ax.invert_yaxis()
    ax.set_xlim(-0.85, 2.55)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend placed cleanly in the middle corridor
    legend_handles = [
        Line2D([0], [0], color=C_A, ls="-", lw=1.5, marker="o", markersize=4.5,
               label="Track A (RNAi-validated benchmark)"),
        Line2D([0], [0], color=C_B, ls="--", lw=1.5, marker="o", markersize=4.5,
               label="Track B (Novel candidate)"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(0.28, 0.68),
              frameon=False, fontsize=6.2)

    fig.suptitle("Candidate Rank Trajectories Across Dirichlet Prior Weighting",
                 fontsize=8.0, y=0.98)
    fig.subplots_adjust(left=0.22, right=0.88, top=0.90, bottom=0.10)
    save(fig, "15_method_bumpchart")


if __name__ == "__main__":
    build()
