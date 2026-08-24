"""Uniform Dirichlet — fixed score vs uniform median score (all 99)."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np
from scipy.stats import spearmanr

def build():
    unif = load_unif99()
    fig, ax = plt.subplots(figsize=(6, 5))
    x = unif["fixed_weight_score"].values
    y = unif["uniform_median_score"].values
    mask = ~(np.isnan(x)|np.isnan(y))
    x, y = x[mask], y[mask]
    rho, p = spearmanr(x, y)
    ax.scatter(x, y, s=12, c=C_UNIFORM, alpha=0.5, edgecolors="white", lw=0.3)
    lo, hi = min(x.min(),y.min())-0.02, max(x.max(),y.max())+0.02
    ax.plot([lo,hi],[lo,hi],"--",color="#999",lw=0.8,label="y = x")
    ax.text(0.05, 0.95, f"Spearman rho = {rho:.3f}\np = {p:.1e}", transform=ax.transAxes,
            fontsize=8, va="top", bbox=dict(boxstyle="round",fc="white",alpha=0.8))
    ax.set_xlabel("Fixed-weight score"); ax.set_ylabel("Uniform Dirichlet median score")
    ax.set_title("Fixed vs uniform Dirichlet score (all 99)", fontweight="bold", pad=8)
    ax.legend(frameon=False, fontsize=7); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "13_uniform_scatter_99")

if __name__=="__main__": build()
