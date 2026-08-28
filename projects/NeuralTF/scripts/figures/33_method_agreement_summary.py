"""Method agreement summary — consensus with significance annotations."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
from matplotlib.patches import Patch

def build():
    try:
        data_path = RES / "overlap_significance.json"
        with open(data_path) as f:
            data = json.load(f)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

        # Panel A: Upset-style overlap counts
        overlaps = data.get("overlaps", {})
        labels = list(overlaps.keys())
        counts = [overlaps[k]["count"] for k in labels]
        pvals = [overlaps[k].get("p_value", None) for k in labels]

        y = np.arange(len(labels))
        colors = [C_A if c > 0 else C_NEURAL for c in counts]
        bars = ax1.barh(y, counts, color=colors, edgecolor="white", lw=0.5, height=0.6)

        for i, (c, p) in enumerate(zip(counts, pvals)):
            sig = ""
            if p is not None:
                sig = f" (p={p:.2e})" if p < 0.05 else " (ns)"
            ax1.text(c + 0.1, i, f"{c}{sig}", va="center", fontsize=7, fontweight="bold")

        ax1.set_yticks(y)
        ax1.set_yticklabels(labels, fontsize=8)
        ax1.set_xlabel("Number of candidates in overlap")
        ax1.set_title("Method overlap counts\nwith significance", fontweight="bold")
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)

        # Panel B: Agreement heatmap (pairwise)
        methods = ["fixed", "centered", "uniform"]
        pairwise = data.get("pairwise", {})
        matrix = np.zeros((3, 3))
        for i, m1 in enumerate(methods):
            for j, m2 in enumerate(methods):
                if m1 == m2:
                    matrix[i, j] = 1.0
                else:
                    key = f"{m1}_vs_{m2}"
                    matrix[i, j] = pairwise.get(key, {}).get("jaccard", 0)

        im = ax2.imshow(matrix, cmap="YlOrRd", vmin=0, vmax=1, aspect="auto")
        ax2.set_xticks(range(3))
        ax2.set_xticklabels(["Fixed", "Centered", "Uniform"], fontsize=8)
        ax2.set_yticks(range(3))
        ax2.set_yticklabels(["Fixed", "Centered", "Uniform"], fontsize=8)

        for i in range(3):
            for j in range(3):
                ax2.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center",
                         fontsize=9, fontweight="bold",
                         color="white" if matrix[i, j] > 0.6 else "black")

        fig.colorbar(im, ax=ax2, shrink=0.8, label="Jaccard similarity")
        ax2.set_title("Pairwise method agreement\n(Jaccard index)", fontweight="bold")

        fig.tight_layout()
        save(fig, "33_method_agreement_summary")
    except FileNotFoundError as e:
        print(f"  [SKIP] {__file__}: {e}")
        return

if __name__ == "__main__":
    build()
