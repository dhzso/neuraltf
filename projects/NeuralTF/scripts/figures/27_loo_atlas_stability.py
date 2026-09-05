"""Leave-one-STREAM-out stability heatmap for top-10 candidates.

The atlas contributions are collapsed into single evidence streams in
rank.csv, so this ablation removes one STREAM (of the 9), not one atlas.
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    try:
        df = pd.read_csv(RES / "loo_atlas_stability.csv")
        streams = [c for c in df.columns if c not in ("gene_id", "gene_name", "full_rank", "full_score")]
        gene_col = "gene_name" if "gene_name" in df.columns else "gene_id"
        names = df[gene_col].values
        rank_matrix = df[streams].values

        fig, ax = plt.subplots(figsize=(8, 5))
        cmap = plt.cm.RdYlGn_r
        im = ax.imshow(rank_matrix, cmap=cmap, aspect="auto", vmin=1, vmax=10)

        ax.set_xticks(range(len(streams)))
        ax.set_xticklabels([a.replace("_", " ").title() for a in streams], rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names, fontsize=8, fontweight="bold")

        for i in range(len(names)):
            for j in range(len(streams)):
                val = rank_matrix[i, j]
                color = "white" if val > 6 else "black"
                ax.text(j, i, f"{int(val)}", ha="center", va="center", fontsize=7, color=color, fontweight="bold")

        cbar = fig.colorbar(im, ax=ax, shrink=0.8, label="Rank when stream is excluded")
        ax.set_title("Top-10 stability across leave-one-stream-out ablations\n"
                     "Lower rank (greener) = more stable despite stream removal",
                     fontweight="bold", pad=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        save(fig, "27_loo_atlas_stability")
    except FileNotFoundError as e:
        print(f"  [SKIP] {__file__}: {e}")
        return

if __name__ == "__main__":
    build()
