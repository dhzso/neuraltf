"""Pairwise Spearman correlation matrix across all 11 evidence streams for the transcriptome."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib.colors as mcolors

def build():
    all_cand = load_all()
    mat = all_cand[STREAM_COLS].fillna(0).values
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

    fig, ax = plt.subplots(figsize=(W_15COL, 4.0))
    im = ax.imshow(corr, cmap="RdBu_r", norm=mcolors.TwoSlopeNorm(vmin=-0.4, vcenter=0.0, vmax=1.0), aspect="equal")
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=6.5)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=6.5)

    # Annotate cells with adaptive contrast text color
    for i in range(n):
        for j in range(n):
            v = corr[i, j]
            if i == j:
                ax.text(j, i, "1.00", ha="center", va="center", fontsize=5.8, color="white")
            else:
                sig = "***" if pvals[i, j] < 0.001 else ("**" if pvals[i, j] < 0.01 else ("*" if pvals[i, j] < 0.05 else ""))
                tc = "white" if (v > 0.48 or v < -0.3) else "#222222"
                ax.text(j, i, f"{v:.2f}{sig}", ha="center", va="center", fontsize=5.5, color=tc)

    title_block(
        fig,
        f"Pairwise Stream Correlation Across Transcriptome (Spearman \u03c1, N = {len(all_cand):,})",
        "*** $P$ < 0.001; ** $P$ < 0.01; * $P$ < 0.05",
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Spearman \u03c1", fontsize=6.8)
    cbar.ax.tick_params(labelsize=6.0)
    ax.spines[:].set_visible(False)
    fig.tight_layout()
    fig.subplots_adjust(top=0.86, bottom=0.16)
    save(fig, "20_stream_correlation")

if __name__=="__main__": build()

