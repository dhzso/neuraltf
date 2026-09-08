"""3-method comparison — rank correlation heatmap."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from scipy.stats import spearmanr

def build():
    s1 = pd.read_csv(RES / "supplementary_table_S1_method_comparison.csv")
    
    methods = ["fixed_composite", "centered_composite", "uniform_composite"]
    labels = ["Fixed-weight", "Centered\nDirichlet", "Uniform\nDirichlet"]
    n = len(methods)
    corr = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            if i == j:
                corr[i, j] = 1.0
            elif i < j:
                rho, _ = spearmanr(s1[methods[i]], s1[methods[j]])
                corr[i, j] = rho
                corr[j, i] = rho

    fig, ax = plt.subplots(figsize=(W_1COL, 3.0))
    cmap = plt.cm.Blues
    im = ax.imshow(corr, cmap=cmap, vmin=0.90, vmax=1.0, aspect="equal")
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=7)

    for i in range(n):
        for j in range(n):
            val = corr[i, j]
            tc = "white" if val > 0.96 else "#222222"
            ax.text(j, i, f"{val:.3f}", ha="center", va="center", fontsize=7.5,
                    color=tc, fontweight="bold")
                    
    ax.set_title("Rank correlation (Spearman $r_s$)",
                 fontweight="bold", fontsize=8.5, pad=6)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.06)
    cbar.set_label("Spearman $r_s$", fontsize=7.5)
    cbar.ax.tick_params(labelsize=6.5)
    
    ax.spines[:].set_visible(False)
    fig.tight_layout()
    save(fig, "17_method_rank_correlation")

if __name__=="__main__": build()

