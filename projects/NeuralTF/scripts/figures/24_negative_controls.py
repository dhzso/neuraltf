"""Negative controls — violin/box plot: neural TFs vs non-TFs vs random."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import json

def build():
    try:
        data_path = RES / "negative_control_stats.json"
        with open(data_path) as f:
            data = json.load(f)

        fig, ax = plt.subplots(figsize=(6, 5))

        groups = ["neural_tfs", "non_tfs", "random"]
        # Labels match the 2026-09-06 stats script: LABEL-FREE score
        # (rnai/neural_*/perez_lineage excluded — the rnai stream is a
        # perfect copy of the group label), availability-matched, and the
        # TF controls are additionally neural_enriched-free.
        labels = ["Neural TFs\n(RNAi-validated)",
                  "No TF class\n(matched controls)",
                  "TF-classified,\nmatched, non-neural"]
        colors = [C_A, C_B, C_NEURAL]
        positions = [1, 2, 3]

        for pos, group, label, color in zip(positions, groups, labels, colors):
            scores = np.array(data[group]) if data.get(group) else np.array([])
            if len(scores) == 0:
                # empty group (e.g. no eligible controls) — never plot the
                # old [0.0] sentinel as a phantom gene
                continue
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
        ax.set_ylabel("Label-free score (rnai/neural/lineage streams excluded)")
        ax.set_title("Neural TFs score higher than availability-matched controls\n(label-free score; controls also neural-enrichment-free)",
                     fontweight="bold", pad=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # 2026-09-06: show BOTH contrasts' p-values — the old figure showed
        # only the weaker non-TF comparison and hid the stricter TF one.
        pvals = []
        if "neural_vs_random_non_tf" in data:
            pvals.append(f"vs non-TF: p={data['neural_vs_random_non_tf']['p_value']:.2e}")
        if "neural_vs_non_neural_tf" in data:
            pvals.append(f"vs non-neural TF: p={data['neural_vs_non_neural_tf']['p_value']:.2e}")
        if pvals:
            ax.annotate("\n".join(pvals),
                        xy=(2.0, max(data["neural_tfs"]) * 1.05) if data.get("neural_tfs") else (2.0, 1.0),
                        ha="center", fontsize=8, fontweight="bold", color=C_HL)

        fig.tight_layout()
        save(fig, "24_negative_controls")
    except FileNotFoundError as e:
        print(f"  [SKIP] {__file__}: {e}")
        return

if __name__ == "__main__":
    build()
