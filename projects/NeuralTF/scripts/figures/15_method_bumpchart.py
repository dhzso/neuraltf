"""Candidate rank stability trajectories across Dirichlet prior weighting schemes.

Dual-track bump chart tracking within-track rank conservation from Fixed weights to
Centered Dirichlet (k=40) and Uniform Dirichlet (alpha=1) for Track A and Track B candidates.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def get_track_ranks(df, score_col="composite_score"):
    if "gene_id_v6" in df.columns and "gene_id" not in df.columns:
        df = df.rename(columns={"gene_id_v6": "gene_id"})
    df_sorted = df.sort_values(score_col, ascending=False).reset_index(drop=True)
    df_a = df_sorted[df_sorted["proof_status"] == "tested"].reset_index(drop=True)
    df_a["rank_a"] = df_a.index + 1
    df_b = df_sorted[df_sorted["proof_status"] == "not_tested"].reset_index(drop=True)
    df_b["rank_b"] = df_b.index + 1
    return df_a.set_index("gene_id")["rank_a"], df_b.set_index("gene_id")["rank_b"]


def build():
    top10_fixed = load_top10()
    top10_cen = load_centered()
    top10_uni = load_uniform()
    cf = load_centered_full()
    uf = load_uniform_full()
    
    fixed_path = RES / "fixed_full_rank.csv"
    if fixed_path.exists():
        fixed_full = pd.read_csv(fixed_path)
        if "gene_id_v6" in fixed_full.columns and "gene_id" not in fixed_full.columns:
            fixed_full = fixed_full.rename(columns={"gene_id_v6": "gene_id"})
    else:
        neural = load_neural()
        neural_bonus = cf.set_index("gene_id")["bonus_total"].reindex(neural["gene_id"]).fillna(0.0).values
        neural["composite_score"] = neural["integrated_score"] + neural_bonus
        fixed_full = neural

    # Rank lookups across the 3 methods using within-track composite scores
    a_fix, b_fix = get_track_ranks(fixed_full, score_col="composite_score")
    a_cen, b_cen = get_track_ranks(cf, score_col="composite_score")
    a_uni, b_uni = get_track_ranks(uf, score_col="composite_score")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.8, 4.0), dpi=500)
    methods = ["Fixed\nWeight", "Centered\nDirichlet (k=40)", "Uniform\nDirichlet (\u03b1=1)"]
    x_pos = [0, 1, 2]

    palette_a = ["#1B365D", "#2B4C6F", "#4A7C59", "#7D5A7D", "#C08A3E", "#65799B"]
    palette_b = ["#B04A3E", "#D9822B", "#5C82A6", "#8C564B", "#2CA02C", "#9467BD", "#E377C2"]

    # ------------------ PANEL A: Tested (RNAi-screened benchmark) ------------------
    panel_tag(ax1, "a", x=-0.14, y=1.05)
    ax1.set_title("Tested (RNAi-screened benchmark) regulators\u2020", fontsize=7.5, pad=8, fontweight="bold")

    # Union of Track A candidates across all 3 methods
    gid_col_f = "gene_id" if "gene_id" in top10_fixed.columns else "gene_id_v6"
    gid_col_c = "gene_id" if "gene_id" in top10_cen.columns else "gene_id_v6"
    gid_col_u = "gene_id" if "gene_id" in top10_uni.columns else "gene_id_v6"

    genes_a = list(dict.fromkeys(
        list(top10_fixed[top10_fixed["track"] == "A"][gid_col_f]) +
        list(top10_cen[top10_cen["track"] == "A"][gid_col_c]) +
        list(top10_uni[top10_uni["track"] == "A"][gid_col_u])
    ))

    for idx, gid in enumerate(genes_a):
        row = cf[cf["gene_id"] == gid].iloc[0] if len(cf[cf["gene_id"] == gid]) else fixed_full[fixed_full["gene_id"] == gid].iloc[0]
        nm = clean_gene_symbol(row.get("gene_name", ""), gid)
        r0 = a_fix.get(gid, 8)
        r1 = a_cen.get(gid, 8)
        r2 = a_uni.get(gid, 8)
        c = palette_a[idx % len(palette_a)]
        ls = "-" if r0 <= 5 and r2 <= 5 else "--"

        ax1.plot(x_pos, [r0, r1, r2], marker="o", color=c, ls=ls, lw=1.6, markersize=5.5, zorder=4)
        ax1.text(-0.08, r0, f"#{r0} {nm}", ha="right", va="center", fontsize=6.2, fontweight="bold", color=c)
        ax1.text(2.08, r2, f"#{r2} {nm}", ha="left", va="center", fontsize=6.2, fontweight="bold", color=c)

    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(methods, fontsize=6.8)
    ax1.set_ylabel("Prioritization Rank (Tested)", fontsize=7.0, fontweight="bold")
    ax1.set_ylim(0.5, 5.8)
    ax1.invert_yaxis()
    ax1.set_xlim(-0.75, 2.75)
    ax1.tick_params(left=False, labelleft=False)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_visible(False)
    ax1.grid(axis="y", color="#EEEEEE", lw=0.6, ls=":")

    # ------------------ PANEL B: Not tested ------------------
    panel_tag(ax2, "b", x=-0.14, y=1.05)
    ax2.set_title("Not tested (no RNAi record) neural TF candidates", fontsize=7.5, pad=8, fontweight="bold")

    genes_b = list(dict.fromkeys(
        list(top10_fixed[top10_fixed["track"] == "B"][gid_col_f]) +
        list(top10_cen[top10_cen["track"] == "B"][gid_col_c]) +
        list(top10_uni[top10_uni["track"] == "B"][gid_col_u])
    ))

    for idx, gid in enumerate(genes_b):
        row = cf[cf["gene_id"] == gid].iloc[0] if len(cf[cf["gene_id"] == gid]) else fixed_full[fixed_full["gene_id"] == gid].iloc[0]
        nm = clean_gene_symbol(row.get("gene_name", ""), gid)
        r0 = b_fix.get(gid, 8)
        r1 = b_cen.get(gid, 8)
        r2 = b_uni.get(gid, 8)
        c = palette_b[idx % len(palette_b)]
        ls = "-" if r0 <= 5 and r2 <= 5 else "--"

        ax2.plot(x_pos, [r0, r1, r2], marker="s", color=c, ls=ls, lw=1.6, markersize=5.0, zorder=4)
        suffix = " (entrant)" if r0 > 5 and r2 <= 5 else (" (displaced)" if r0 <= 5 and r2 > 5 else "")
        ax2.text(-0.08, r0, f"#{r0} {nm}", ha="right", va="center", fontsize=6.2, fontweight="bold", color=c)
        ax2.text(2.08, r2, f"#{r2} {nm}{suffix}", ha="left", va="center", fontsize=6.2, fontweight="bold", color=c)

    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(methods, fontsize=6.8)
    ax2.set_ylabel("Prioritization Rank (Not tested)", fontsize=7.0, fontweight="bold")
    ax2.set_ylim(0.5, 6.8)
    ax2.invert_yaxis()
    ax2.set_xlim(-0.75, 2.95)
    ax2.tick_params(left=False, labelleft=False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.grid(axis="y", color="#EEEEEE", lw=0.6, ls=":")

    fig.suptitle("Prioritization Rank Trajectories Across Dirichlet Prior Weighting Schemes",
                 fontsize=8.5, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.10, right=0.88, top=0.84, bottom=0.14, wspace=0.42)
    save(fig, "15_method_bumpchart")


if __name__ == "__main__":
    build()
