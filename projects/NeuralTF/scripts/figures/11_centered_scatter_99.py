"""Centered Dirichlet — fixed-weight score vs centered median (top 10 only).
NOTE: Centered Dirichlet script only saved top-10 results, not all 99."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from scipy.stats import spearmanr

def build():
    c = load_centered()
    fig, ax = plt.subplots(figsize=(6, 5))
    x = c["dirichlet_median_score"].values
    y = c["composite_score"].values
    tracks = c["track"].values
    names = [label(load_neural(), g) for g in c["gene_id"]]
    for i in range(len(x)):
        tc = C_A if tracks[i]=="A" else C_B
        ax.scatter(x[i], y[i], s=50, c=tc, edgecolors="white", lw=0.5, zorder=5)
        ax.annotate(names[i], (x[i], y[i]), fontsize=7, ha="center", va="bottom",
                    xytext=(0,5), textcoords="offset points", color=tc, fontweight="bold")
    lo = min(x.min(), y.min()) - 0.02
    hi = max(x.max(), y.max()) + 0.02
    ax.plot([lo,hi],[lo,hi],"--",color="#999",lw=0.8,label="y = x")
    rho, p = spearmanr(x, y)
    ax.text(0.05, 0.95, f"n = {len(x)} (top-10 only)\nSpearman rho = {rho:.3f}\np = {p:.1e}",
            transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(boxstyle="round",fc="white",alpha=0.8))
    ax.set_xlabel("Centered Dirichlet median score"); ax.set_ylabel("Composite score (median + bonuses)")
    ax.set_title("Centered Dirichlet — median vs composite (top 10)\n"
                 "Note: full 99-candidate centered data not saved by pipeline",
                 fontweight="bold", pad=8, fontsize=9)
    ax.legend(frameon=False, fontsize=7); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "11_centered_scatter_99")

if __name__=="__main__": build()
