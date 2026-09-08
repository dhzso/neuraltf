"""Uniform Dirichlet — integrated score vs uniform median (all TF candidates)."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from scipy.stats import spearmanr

def build():
    centered = load_centered_full()
    uniform = load_uniform_full()
    all_df = load_all()

    # Merge on gene_id
    df = all_df[["gene_id", "integrated_score", "proof_status"]].merge(
        centered[["gene_id", "dirichlet_median_score"]], on="gene_id", how="inner"
    ).merge(
        uniform[["gene_id", "uniform_median_score"]], on="gene_id", how="inner"
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_2COL, 3.4), sharey=True)

    x = df["integrated_score"].values
    proof = df["proof_status"].fillna("").values
    colors = [C_A if "validated" in str(p).lower() else (C_B if "novel" in str(p).lower() else "#B0BEC5") for p in proof]

    panels = [
        (ax1, df["dirichlet_median_score"].values, "Centered Dirichlet (k=40)", "a"),
        (ax2, df["uniform_median_score"].values, "Uniform Dirichlet (\u03b1=1)", "b")
    ]

    for ax, y, title, tag in panels:
        mask = ~(np.isnan(x) | np.isnan(y))
        x_m, y_m = x[mask], y[mask]
        
        # Background candidates
        bg_mask = np.array([c == "#B0BEC5" for c in colors])[mask]
        ax.scatter(x_m[bg_mask], y_m[bg_mask], s=8, color="#B0BEC5", alpha=0.35, edgecolors="none")
        
        # Highlight candidates (Track A & B)
        hl_mask = ~bg_mask
        ax.scatter(x_m[hl_mask], y_m[hl_mask], s=22, c=np.array(colors)[mask][hl_mask],
                   alpha=0.9, edgecolors="white", lw=0.4, zorder=5)

        lo, hi = -0.02, 1.05
        ax.plot([lo, hi], [lo, hi], "--", color="#555555", lw=0.8, label="y = x (identity)")
        
        rho, p = spearmanr(x_m, y_m)
        p_str = "p < 10^{-30}" if p < 1e-30 else f"p = {p:.1e}"
        ax.text(0.06, 0.92, f"$r_s = {rho:.3f}$\n${p_str}$\n(n = {len(x_m):,})",
                transform=ax.transAxes, fontsize=7.5, va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="#FAFAFA", ec="#CCCCCC", lw=0.5))

        ax.set_xlabel("Fixed-weight integrated score", fontsize=8)
        ax.set_title(title, fontweight="bold", fontsize=8.5, pad=6)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        panel_tag(ax, tag)

    ax1.set_ylabel("Dirichlet median score", fontsize=8)

    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_A, markersize=6, label="Track A (RNAi-validated)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_B, markersize=6, label="Track B (novel candidate)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#B0BEC5", markersize=5, label="Other candidate"),
        Line2D([0], [0], color="#555555", ls="--", lw=0.8, label="Identity line")
    ]
    ax2.legend(handles=legend_handles, frameon=False, fontsize=6.8, loc="lower right")

    fig.suptitle("Prioritization concordance: fixed-weight vs Bayesian Dirichlet scoring",
                 fontweight="bold", fontsize=9, y=1.01)
    fig.tight_layout()
    save(fig, "13_uniform_scatter_all")

if __name__=="__main__": build()

