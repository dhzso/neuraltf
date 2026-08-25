"""Centered Dirichlet — fixed-weight score vs centered median (all 99 candidates)."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np
from scipy.stats import spearmanr

def build():
    full = load_centered99()
    neural = load_neural()
    fig, ax = plt.subplots(figsize=(6, 5))
    x = full["fixed_weight_score"].values
    y = full["dirichlet_median_score"].values
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    proof = full.loc[mask, "proof_status"].values if "proof_status" in full.columns else [""]*len(x)
    colors = [C_A if "validated" in str(p).lower() else C_B if "novel" in str(p).lower() else C_NEURAL for p in proof]
    ax.scatter(x, y, s=15, c=colors, alpha=0.5, edgecolors="white", lw=0.3)
    lo, hi = min(x.min(),y.min())-0.02, max(x.max(),y.max())+0.02
    ax.plot([lo,hi],[lo,hi],"--",color="#999",lw=0.8,label="y = x")
    rho, p = spearmanr(x, y)
    ax.text(0.05, 0.95, f"n = {len(x)} candidates\nSpearman rho = {rho:.3f}\np = {p:.1e}",
            transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(boxstyle="round",fc="white",alpha=0.8))
    ax.set_xlabel("Fixed-weight integrated score"); ax.set_ylabel("Centered Dirichlet median score (k=40)")
    ax.set_title(f"Fixed-weight and centered Dirichlet scores are strongly correlated (rho = {rho:.2f})",
                 fontweight="bold", pad=8)
    ax.legend(frameon=False, fontsize=7); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "11_centered_scatter_99")

if __name__=="__main__": build()
