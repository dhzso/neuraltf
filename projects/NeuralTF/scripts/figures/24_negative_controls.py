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

    fig, ax = plt.subplots(figsize=(5.4, 4.0), dpi=500)

    groups = ["neural_tfs", "non_tfs", "random"]
    labels = [
        "Neural TFs\n(RNAi+, n=67)",
        "Matched Non-TFs\n(n=100)",
        "Matched Non-neural TFs\n(n=100)",
    ]
    colors = [C_A, "#788896", C_B]
    positions = [1, 2, 3]

    np.random.seed(42)
    final_labels = []

    for pos, group, label_str, color in zip(positions, groups, labels, colors):
        scores = np.array(data[group]) if data.get(group) else np.array([])
        if len(scores) == 0:
            continue

        # Violin density body
        parts = ax.violinplot(scores, positions=[pos], showmeans=False, showmedians=False, widths=0.65)
        for pc in parts["bodies"]:
            pc.set_facecolor(color)
            pc.set_edgecolor("none")
            pc.set_alpha(0.25)
        for key in ("cbars", "cmins", "cmaxes"):
            if key in parts:
                parts[key].set_visible(False)

        # Jittered data points overlay
        jitter = np.random.normal(0, 0.04, size=len(scores))
        ax.scatter(pos + jitter, scores, s=12, color=color, alpha=0.35, edgecolors="none", zorder=3)

        # Boxplot overlay
        bp = ax.boxplot(scores, positions=[pos], widths=0.22, patch_artist=True,
                        showfliers=False, zorder=4,
                        boxprops=dict(facecolor=color, alpha=0.85, edgecolor="none"),
                        medianprops=dict(color="white", lw=1.5),
                        whiskerprops=dict(color="#444444", lw=0.8),
                        capprops=dict(color="#444444", lw=0.8))

        median = np.median(scores)
        final_labels.append(f"{label_str}\n(med = {median:.2f})")

    ax.set_xticks(positions)
    ax.set_xticklabels(final_labels, fontsize=6.8)
    ax.set_ylabel("Label-Free Evidence Score (0–1)", fontsize=7.0, fontweight="bold")
    ax.set_title("Biological Specificity: Score Distribution vs Empirical Negative Controls",
                 fontsize=8.0, pad=10, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0.35, 3.65)

    def _format_p(p):
        if p is None:
            return ""
        if p < 1e-15:
            return r"P < 10^{-15}"
        base, exp = f"{p:.1e}".split("e")
        return rf"P = {base} \times 10^{{{int(exp)}}}"

    # Statistical contrasts from JSON
    s1 = data.get("neural_vs_random_non_tf", {})
    s2 = data.get("neural_vs_non_neural_tf", {})
    u1 = s1.get("mann_whitney_u", None)
    u2 = s2.get("mann_whitney_u", None)
    p1 = s1.get("p_value", None)
    d1 = s1.get("cohens_d", None)
    p2 = s2.get("p_value", None)
    d2 = s2.get("cohens_d", None)

    y_bar1 = 1.05
    h = 0.025
    if p1 is not None:
        ax.plot([1, 1, 2, 2], [y_bar1, y_bar1 + h, y_bar1 + h, y_bar1], color="#333333", lw=0.7)
        p1_str = _format_p(p1)
        d1_str = f", $d = {d1:.2f}$" if d1 else ""
        ax.text(1.5, y_bar1 + h + 0.015, f"${p1_str}${d1_str}", ha="center", va="bottom", fontsize=6.2, color="#222222")

    y_bar2 = 1.18
    if p2 is not None:
        ax.plot([1, 1, 3, 3], [y_bar2, y_bar2 + h, y_bar2 + h, y_bar2], color="#333333", lw=0.7)
        p2_str = _format_p(p2)
        d2_str = f", $d = {d2:.2f}$" if d2 else ""
        ax.text(2.0, y_bar2 + h + 0.015, f"${p2_str}${d2_str}", ha="center", va="bottom", fontsize=6.2, color="#222222")

    # Footnote note explaining test, effect size, circularity control & matching
    u1_str = f"{int(round(u1)):,}" if u1 is not None else "5,577"
    u2_str = f"{int(round(u2)):,}" if u2 is not None else "5,236"
    ax.text(0.5, -0.16,
            f"Two-sided Mann–Whitney U test ($U_1 = {u1_str}, U_2 = {u2_str}$) and Cohen's d effect size.\n"
            "Label-free score excludes RNAi, neural enrichment & neural lineage to eliminate circularity.\n"
            "Control cohorts matched on number of available evidence streams.",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.8, color="#555555", style="italic")

    ax.set_ylim(-0.02, 1.34)
    fig.subplots_adjust(left=0.12, right=0.96, top=0.90, bottom=0.22)
    save(fig, "24_negative_controls")

if __name__ == "__main__":
    build()
