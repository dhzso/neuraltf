"""Uniform Dirichlet robustness — integrated score vs uninformative Dirichlet median (all candidates).

Single-panel scatter plot evaluating global score stability across all planarian TFs under
1,000 Dirichlet weight draws from a uniform prior (alpha = 1).
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
from matplotlib.lines import Line2D


def build():
    uniform = load_uniform_full()
    all_df = load_all()

    # Merge on gene_id
    df = all_df[["gene_id", "integrated_score", "proof_status"]].merge(
        uniform[["gene_id", "uniform_median_score"]], on="gene_id", how="inner"
    )

    fig, ax = plt.subplots(figsize=(W_15COL, 3.8), dpi=500)

    x = df["integrated_score"].values
    y = df["uniform_median_score"].values
    proof = df["proof_status"].fillna("").values

    mask = ~(np.isnan(x) | np.isnan(y))
    x_m, y_m = x[mask], y[mask]
    proof_m = proof[mask]

    neural = load_neural()
    track_a_ids = set(neural[neural["proof_status"] == "known_rnai_validated"]["gene_id"])
    track_b_ids = set(neural[neural["proof_status"] == "novel_candidate"]["gene_id"])

    # Partition into Background, Track A, Track B
    is_val = np.array([gid in track_a_ids for gid in df["gene_id"]])[mask]
    is_nov = np.array([gid in track_b_ids for gid in df["gene_id"]])[mask]
    is_bg = ~(is_val | is_nov)

    # 1. Background candidates (unfiltered transcriptome)
    ax.scatter(x_m[is_bg], y_m[is_bg], s=8, color="#C8CED6", alpha=0.3,
               edgecolors="none", label=f"Transcriptome-wide (n={np.sum(is_bg):,})")

    # 2. Track A (RNAi-validated benchmark)
    ax.scatter(x_m[is_val], y_m[is_val], s=28, color=C_A, alpha=0.92,
               edgecolors="white", lw=0.5, zorder=5, label=f"Track A: Validated (n={np.sum(is_val)})")

    # 3. Track B (Novel candidates)
    ax.scatter(x_m[is_nov], y_m[is_nov], s=28, color=C_B, alpha=0.92,
               edgecolors="white", lw=0.5, zorder=6, label=f"Track B: Novel (n={np.sum(is_nov)})")

    # Identity reference line
    lo, hi = -0.02, 1.05
    ax.plot([lo, hi], [lo, hi], "--", color="#666666", lw=0.8, label="y = x (identity)")

    # Spearman correlation
    rho, p = spearmanr(x_m, y_m)
    p_str = "P < 10^{-300}" if p < 1e-300 else f"P = {p:.1e}"
    ax.text(0.05, 0.93, f"Spearman $r_s = {rho:.3f}$\n${p_str}$\n$N = {len(x_m):,}$",
            transform=ax.transAxes, fontsize=7.5, va="top", color="#222222",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#DDDDDD", lw=0.5))

    ax.set_xlabel("Fixed Weight Integrated Score", fontsize=8, fontweight="bold")
    ax.set_ylabel("Uniform Dirichlet Median Score (1,000 Draws)", fontsize=8, fontweight="bold")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend cleanly positioned at lower right
    ax.legend(loc="lower right", frameon=False, fontsize=7)

    fig.suptitle("Score Invariance Under Uninformative Dirichlet Prior Weighting",
                 fontweight="bold", fontsize=8.5, y=0.98)
    fig.subplots_adjust(left=0.14, right=0.96, top=0.90, bottom=0.14)
    save(fig, "13_uniform_scatter_all")


if __name__ == "__main__":
    build()
