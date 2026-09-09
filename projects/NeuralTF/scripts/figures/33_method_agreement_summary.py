"""Method agreement summary — consensus with significance annotations.

Panel a: Method overlap counts with hypergeometric test significance.
Panel b: Pairwise method agreement (Jaccard similarity matrix).
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import json

def build():
    data_path = RES / "overlap_significance.json"
    if not data_path.exists():
        raise FileNotFoundError(
            f"{data_path} missing — run scripts/stats/overlap_significance.py first"
        )
    with open(data_path) as f:
        data = json.load(f)

    # 3x3 Pairwise Jaccard similarity matrix across prioritization formulations
    methods = ["fixed", "centered", "uniform"]
    display_methods = [
        "Fixed weights",
        "Dirichlet Centered\n($k = 40$)",
        "Dirichlet Uniform\n($\\alpha = 1$)",
    ]
    pairwise = data.get("pairwise", {})
    matrix = np.zeros((3, 3))
    counts = np.zeros((3, 3), dtype=int)

    for i, m1 in enumerate(methods):
        for j, m2 in enumerate(methods):
            if i == j:
                matrix[i, j] = 1.0
                counts[i, j] = 10
            else:
                key = f"{m1}_vs_{m2}"
                info = pairwise.get(key, {})
                matrix[i, j] = info.get("jaccard", 0.0)
                counts[i, j] = info.get("overlap_count", 0)

    fig, ax = plt.subplots(figsize=(3.8, 3.3))

    # Muted blue sequential colormap
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "custom_blues", ["#F5F8FA", "#D4E2EE", "#7AA2C0", "#2B4C6F"], N=256
    )

    im = ax.imshow(matrix, cmap=cmap, vmin=0.0, vmax=1.0, aspect="auto")

    ax.set_xticks(range(3))
    ax.set_xticklabels(display_methods, fontsize=6.8, rotation=25, ha="right")
    ax.set_yticks(range(3))
    ax.set_yticklabels(display_methods, fontsize=6.8)

    # Cell annotations with both Jaccard index and candidate count
    for i in range(3):
        for j in range(3):
            val = matrix[i, j]
            cnt = counts[i, j]
            text_color = "white" if val > 0.60 else "#222222"
            if i == j:
                cell_text = f"$J = 1.00$\n(10/10)"
            else:
                cell_text = f"$J = {val:.2f}$\n({cnt}/10 shared)"
            ax.text(
                j,
                i,
                cell_text,
                ha="center",
                va="center",
                fontsize=7.0,
                color=text_color,
            )

    ax.set_title("Prioritization method agreement (top 10 candidates)", fontsize=8, pad=8)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.82, pad=0.04)
    cbar.set_label("Jaccard similarity index", fontsize=7)
    cbar.ax.tick_params(labelsize=6.5)

    # Subtle grid lines between cells
    ax.set_xticks(np.arange(-0.5, 3, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 3, 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    fig.tight_layout()
    save(fig, "33_method_agreement_summary")

if __name__ == "__main__":
    build()

