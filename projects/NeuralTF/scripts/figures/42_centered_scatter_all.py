"""Centered Dirichlet robustness — integrated score vs centered-Dirichlet median (all candidates).

Single-panel scatter plot evaluating global score stability across all planarian TFs under
1,000 Dirichlet weight draws from the concentrated prior (alpha = 40 * W; k = 40).
Companion panel to fig 13 (uniform prior): together they bracket the plausible prior space.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def build():
    centered = load_centered_full()
    all_df = load_all()

    df = all_df[["gene_id", "integrated_score", "proof_status"]].merge(
        centered[["gene_id", "dirichlet_median_score"]], on="gene_id", how="inner"
    )

    fig, ax = plt.subplots(figsize=(W_15COL, 3.8), dpi=500)

    x = df["integrated_score"].values
    y = df["dirichlet_median_score"].values

    mask = ~(np.isnan(x) | np.isnan(y))
    x_m, y_m = x[mask], y[mask]

    neural = load_neural()
    track_a_ids = set(neural[neural["proof_status"] == "tested"]["gene_id"])
    track_b_ids = set(neural[neural["proof_status"] == "not_tested"]["gene_id"])

    is_val = np.array([gid in track_a_ids for gid in df["gene_id"]])[mask]
    is_nov = np.array([gid in track_b_ids for gid in df["gene_id"]])[mask]
    is_bg = ~(is_val | is_nov)

    ax.scatter(x_m[is_bg], y_m[is_bg], s=8, color="#C8CED6", alpha=0.3,
               edgecolors="none", label=f"Transcriptome-wide (n={np.sum(is_bg):,})")
    ax.scatter(x_m[is_val], y_m[is_val], s=28, color=C_A, alpha=0.92,
               edgecolors="white", lw=0.5, zorder=5, label=f"Track A: RNAi-screened (n={np.sum(is_val)})")
    ax.scatter(x_m[is_nov], y_m[is_nov], s=28, color=C_B, alpha=0.92,
               edgecolors="white", lw=0.5, zorder=6, label=f"Track B: not tested (n={np.sum(is_nov)})")

    lo, hi = -0.02, 1.05
    ax.plot([lo, hi], [lo, hi], "--", color="#666666", lw=0.8, label="y = x (identity)")

    rho, p = spearmanr(x_m, y_m)
    p_str = "P < 10^{-300}" if p < 1e-300 else f"P = {p:.1e}"
    top10 = load_top10()
    t = df[df["gene_id"].isin(set(top10["gene_id"]))]
    dmax = float((t["integrated_score"] - t["dirichlet_median_score"]).abs().max())

    ax.set_xlabel("Fixed weight integrated score", fontsize=7.0)
    ax.set_ylabel("Centered Dirichlet median score (1,000 draws)", fontsize=7.0)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=False, fontsize=6.0)

    title_block(
        fig,
        "Score concordance: Fixed weights vs Centered Dirichlet prior (k = 40)",
        f"Spearman $r_s$ = {rho:.3f}, $N$ = {len(x_m):,}; 1,000 concentrated-prior draws; dashed line = identity $y = x$",
    )
    fig.subplots_adjust(left=0.14, right=0.96, top=0.86, bottom=0.15)
    save(fig, "42_centered_scatter_all")


if __name__ == "__main__":
    build()
