"""Evidence stream correlation matrix — how streams correlate across 249 TFs."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from scipy.stats import spearmanr

def build():
    all249 = load_all()
    mat = all249[STREAM_COLS].fillna(0).values
    n = len(STREAM_COLS)
    labels = [STREAM_L[s] for s in STREAM_COLS]
    colors = [STREAM_C[s] for s in STREAM_COLS]

    # Spearman correlation
    corr = np.zeros((n, n))
    pvals = np.ones((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                corr[i,j] = 1.0
            elif i < j:
                rho, p = spearmanr(mat[:,i], mat[:,j])
                corr[i,j] = rho
                corr[j,i] = rho
                pvals[i,j] = p
                pvals[j,i] = p

    fig, ax = plt.subplots(figsize=(W_15COL, 4.2))
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-0.4, vmax=1.0, aspect="equal")
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=7)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=7)

    # Annotate cells
    for i in range(n):
        for j in range(n):
            v = corr[i, j]
            if i == j:
                ax.text(j, i, "1.00", ha="center", va="center", fontsize=6, color="#222222")
            else:
                sig = "***" if pvals[i, j] < 0.001 else ("**" if pvals[i, j] < 0.01 else ("*" if pvals[i, j] < 0.05 else ""))
                tc = "white" if abs(v) > 0.55 else "#222222"
                ax.text(j, i, f"{v:.2f}{sig}", ha="center", va="center", fontsize=5.8, color=tc)

    ax.set_title("Pairwise stream correlation (Spearman $r_s$)",
                 fontweight="bold", fontsize=8.5, pad=8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Spearman $r_s$", fontsize=7.5)
    cbar.ax.tick_params(labelsize=6.5)
    ax.spines[:].set_visible(False)
    fig.tight_layout()
    save(fig, "20_stream_correlation")

if __name__=="__main__": build()

