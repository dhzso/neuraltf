"""Negative controls — violin/box plot: neural TFs vs non-TFs vs random."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import json

def build():
    data_path = RES / "negative_control_stats.json"
    with open(data_path) as f:
        data = json.load(f)

    fig, ax = plt.subplots(figsize=(6, 5))

    groups = ["neural_tfs", "non_tfs", "random"]
    labels = ["Neural TFs\n(validated)", "Non-TFs\n(catalog)", "Random\n(permutation)"]
    colors = [C_A, C_B, C_NEURAL]
    positions = [1, 2, 3]

    for pos, group, label, color in zip(positions, groups, labels, colors):
        scores = np.array(data[group])
        parts = ax.violinplot(scores, positions=[pos], showmeans=True, showmedians=True)
        for pc in parts["bodies"]:
            pc.set_facecolor(color)
            pc.set_alpha(0.4)
        parts["cmeans"].set_color(color)
        parts["cmedians"].set_color(C_HL)
        parts["cbars"].set_color(color)
        parts["cmins"].set_color(color)
        parts["cmaxes"].set_color(color)

        median = np.median(scores)
        q1, q3 = np.percentile(scores, [25, 75])
        ax.text(pos + 0.25, median, f"median={median:.3f}\nIQR=[{q1:.3f}, {q3:.3f}]",
                fontsize=6, va="center")

    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Integrated score")
    ax.set_title("Neural TFs score significantly higher than negative controls",
                 fontweight="bold", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Significance annotations
    if "p_neural_vs_non" in data:
        p = data["p_neural_vs_non"]
        ax.annotate(f"p = {p:.2e}", xy=(1.5, max(data["neural_tfs"]) * 1.05),
                    ha="center", fontsize=8, fontweight="bold", color=C_HL)

    fig.tight_layout()
    save(fig, "24_negative_controls")

if __name__ == "__main__":
    build()
