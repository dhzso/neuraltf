"""Score distribution across all 9 evidence streams — updated evidence heatmap."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    df = load_all()
    df = df.sort_values("composite_score", ascending=False).head(50)

    streams = STREAM_COLS
    matrix = df[streams].values
    names = df["gene_name"].values if "gene_name" in df.columns else df["gene_id"].values

    fig, ax = plt.subplots(figsize=(10, 8))
    cmap = plt.cm.YlOrRd
    im = ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=0, vmax=1)

    ax.set_xticks(range(len(streams)))
    ax.set_xticklabels([STREAM_L[s] for s in streams], rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=7, fontweight="bold")

    # Annotate cells
    for i in range(len(names)):
        for j in range(len(streams)):
            val = matrix[i, j]
            if pd.notna(val) and val > 0:
                color = "white" if val > 0.7 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5, color=color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, label="Evidence score (0–1)")
    ax.set_title("Top-50 candidates: evidence stream scores across all 9 streams\n"
                 "Darker = stronger evidence in that stream",
                 fontweight="bold", pad=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "31_score_distribution_all9")

if __name__ == "__main__":
    build()
