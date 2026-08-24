"""Centered vs Uniform Dirichlet — median score comparison (all 99 candidates)."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np
from scipy.stats import spearmanr

def build():
    centered = load_centered99()
    uniform = load_unif99()
    neural = load_neural()

    merged = centered[["gene_id","dirichlet_median_score","proof_status"]].merge(
        uniform[["gene_id","uniform_median_score"]], on="gene_id", how="inner")
    fig, ax = plt.subplots(figsize=(6, 5))
    x = merged["dirichlet_median_score"].values
    y = merged["uniform_median_score"].values
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    proof = merged.loc[mask, "proof_status"].values
    colors = [C_A if "validated" in str(p).lower() else C_B if "novel" in str(p).lower() else C_NEURAL for p in proof]
    ax.scatter(x, y, s=15, c=colors, alpha=0.5, edgecolors="white", lw=0.3)
    lo, hi = min(x.min(),y.min())-0.02, max(x.max(),y.max())+0.02
    ax.plot([lo,hi],[lo,hi],"--",color="#999",lw=0.8,label="y = x")
    rho, p = spearmanr(x, y)
    ax.text(0.05, 0.95, f"n = {len(x)} candidates\nSpearman rho = {rho:.3f}\np = {p:.1e}",
            transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(boxstyle="round",fc="white",alpha=0.8))
    ax.set_xlabel("Centered Dirichlet median score (k=40)")
    ax.set_ylabel("Uniform Dirichlet median score (alpha=1)")
    ax.set_title("Centered vs uniform Dirichlet (all 99)", fontweight="bold", pad=8)
    from matplotlib.lines import Line2D
    legend_handles = [Line2D([0],[0], marker="o", color="w", markerfacecolor=C_A, markersize=6, label="RNAi-validated"),
                      Line2D([0],[0], marker="o", color="w", markerfacecolor=C_B, markersize=6, label="Novel"),
                      Line2D([0],[0], marker="o", color="w", markerfacecolor=C_NEURAL, markersize=6, label="Other"),
                      plt.Line2D([0],[0], color="#999", ls="--", lw=0.8, label="y = x")]
    ax.legend(handles=legend_handles, frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "21_centered_vs_uniform_scatter")

if __name__=="__main__": build()
