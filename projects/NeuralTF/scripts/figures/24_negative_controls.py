"""Negative controls — violin/box plot: neural TFs vs non-TFs vs random."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import json

def build():
    data_path = RES / "negative_control_stats.json"
    if not data_path.exists():
        raise FileNotFoundError(f"{data_path} missing")
    with open(data_path) as f:
        data = json.load(f)

    fig, ax = plt.subplots(figsize=(W_15COL, 3.6))

    groups = ["neural_tfs", "non_tfs", "random"]
    labels = ["Neural TFs",
              "Non-TF controls",
              "Non-neural TFs"]
    colors = [C_A, "#788896", C_B]
    positions = [1, 2, 3]

    final_labels = []
    for pos, group, label_str, color in zip(positions, groups, labels, colors):
        scores = np.array(data[group]) if data.get(group) else np.array([])
        if len(scores) == 0:
            continue
        parts = ax.violinplot(scores, positions=[pos], showmeans=False, showmedians=False, widths=0.65)
        for pc in parts["bodies"]:
            pc.set_facecolor(color)
            pc.set_edgecolor("none")
            pc.set_alpha(0.35)
        for key in ("cbars", "cmins", "cmaxes"):
            if key in parts:
                parts[key].set_visible(False)

        # Boxplot overlay
        bp = ax.boxplot(scores, positions=[pos], widths=0.22, patch_artist=True,
                        showfliers=False,
                        boxprops=dict(facecolor=color, alpha=0.85, edgecolor="none"),
                        medianprops=dict(color="white", lw=1.5),
                        whiskerprops=dict(color="#555555", lw=0.8),
                        capprops=dict(color="#555555", lw=0.8))

        median = np.median(scores)
        final_labels.append(f"{label_str}\n(med = {median:.2f})")

    ax.set_xticks(positions)
    ax.set_xticklabels(final_labels, fontsize=6.8)
    ax.set_ylabel("Label-free score", fontsize=7.0)
    ax.set_title("Score distribution vs empirical controls", fontsize=8.0, pad=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0.4, 3.8)

    # Clean p-value brackets
    p1 = data.get("neural_vs_random_non_tf", {}).get("p_value", None)
    p2 = data.get("neural_vs_non_neural_tf", {}).get("p_value", None)
    y_max = max([max(data[g]) for g in groups if data.get(g)]) * 1.05

    if p1 is not None:
        ax.plot([1, 1, 2, 2], [y_max, y_max + 0.04, y_max + 0.04, y_max], color="#444444", lw=0.6)
        p1_str = "p < 10^{-10}" if p1 < 1e-10 else f"p = {p1:.1e}"
        ax.text(1.5, y_max + 0.05, f"${p1_str}$", ha="center", va="bottom", fontsize=6.0, color="#333333")

    if p2 is not None:
        ax.plot([1, 1, 3, 3], [y_max + 0.12, y_max + 0.16, y_max + 0.16, y_max + 0.12], color="#444444", lw=0.6)
        p2_str = "p < 10^{-10}" if p2 < 1e-10 else f"p = {p2:.1e}"
        ax.text(2.0, y_max + 0.17, f"${p2_str}$", ha="center", va="bottom", fontsize=6.0, color="#333333")

    ax.set_ylim(-0.02, y_max + 0.25)
    fig.tight_layout()
    save(fig, "24_negative_controls")

if __name__ == "__main__":
    build()
