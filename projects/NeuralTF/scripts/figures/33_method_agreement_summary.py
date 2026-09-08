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

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_2COL, 2.7),
                                   gridspec_kw={"width_ratios": [1.2, 1]})

    # Panel a: Overlap counts (deduplicate symmetric comparisons)
    overlaps = data.get("overlaps", {})
    comparisons = [
        ("Centered vs Uniform", "centered_vs_uniform"),
        ("Fixed vs Centered", "fixed_vs_centered"),
        ("Fixed vs Uniform", "fixed_vs_uniform"),
        ("All Three Methods", "three_way"),
    ]
    labels = []
    counts = []
    pvals = []
    for display_name, key in comparisons:
        if key in overlaps:
            labels.append(display_name)
            counts.append(overlaps[key]["count"])
            pvals.append(overlaps[key].get("p_value"))

    y = np.arange(len(labels))
    colors = [C_A if c > 5 else C_B for c in counts]
    bars = ax1.barh(y, counts, color=colors, edgecolor="none", height=0.55)

    for i, (c, p) in enumerate(zip(counts, pvals)):
        sig = f" (p={p:.1e})" if (p is not None and p < 0.05) else ""
        ax1.text(c + 0.2, i, f"{c}/10{sig}", va="center", ha="left",
                 fontsize=6.5, fontweight="bold", color="#222222")

    ax1.set_yticks(y)
    ax1.set_yticklabels(labels, fontsize=7)
    ax1.set_xlabel("Candidates in top 10 overlap", fontsize=7.5)
    ax1.set_xlim(0, 14)
    ax1.invert_yaxis()
    panel_tag(ax1, "a")

    # Panel b: Pairwise Jaccard similarity matrix
    methods = ["fixed", "centered", "uniform"]
    display_methods = ["Fixed", "Centered", "Uniform"]
    pairwise = data.get("pairwise", {})
    matrix = np.zeros((3, 3))
    for i, m1 in enumerate(methods):
        for j, m2 in enumerate(methods):
            if i == j:
                matrix[i, j] = 1.0
            else:
                key = f"{m1}_vs_{m2}"
                matrix[i, j] = pairwise.get(key, {}).get("jaccard", 0.0)

    im = ax2.imshow(matrix, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax2.set_xticks(range(3))
    ax2.set_xticklabels(display_methods, fontsize=7)
    ax2.set_yticks(range(3))
    ax2.set_yticklabels(display_methods, fontsize=7)

    for i in range(3):
        for j in range(3):
            val = matrix[i, j]
            ax2.text(j, i, f"{val:.2f}", ha="center", va="center",
                     fontsize=8, fontweight="bold",
                     color="white" if val > 0.6 else "#222222")

    cbar = fig.colorbar(im, ax=ax2, shrink=0.85, pad=0.04)
    cbar.set_label("Jaccard similarity", fontsize=7)
    cbar.ax.tick_params(labelsize=6.5)
    panel_tag(ax2, "b")

    fig.tight_layout()
    save(fig, "33_method_agreement_summary")

if __name__ == "__main__":
    build()

