"""Effect sizes — Cliff's delta + Hedges' g bar panel with honesty labels."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import json

# human-readable labels incl. circularity status
LABELS = {
    "top10_vs_rest": ("Top-10 vs rest (circular)", True),
    "neural_vs_non_neural": ("Neural vs rest (circular)", True),
    "neural_vs_non_neural_honest": ("Neural vs rest (honest)", False),
}

def build():
    try:
        data_path = RES / "effect_sizes.json"
        with open(data_path) as f:
            data = json.load(f)

        comparisons = [k for k in LABELS if k in data]
        if not comparisons:
            raise ValueError("effect_sizes.json carries no known comparisons")

        fig, ax = plt.subplots(figsize=(8, 5))

        y = np.arange(len(comparisons))
        deltas = [data[k]["cliffs_delta"] for k in comparisons]
        gs = [data[k].get("hedges_g", data[k].get("cohens_d", 0.0)) for k in comparisons]
        pvals = [data[k].get("p_value", None) for k in comparisons]

        width = 0.38
        bars_d = ax.barh(y + width/2, deltas, height=width, color=C_A,
                         edgecolor="white", lw=0.5, label="Cliff's delta")
        bars_g = ax.barh(y - width/2, gs, height=width, color=C_B,
                         edgecolor="white", lw=0.5, label="Hedges' g")

        for i, (d, g, p) in enumerate(zip(deltas, gs, pvals)):
            sig = ""
            if p is not None:
                sig = " ***" if p < 0.001 else " **" if p < 0.01 else " *" if p < 0.05 else " ns"
            ax.text(max(d, g) + 0.04, i,
                    f"delta={d:.3f}{sig}", va="center", ha="left",
                    fontsize=7, fontweight="bold")

        ax.set_yticks(y)
        ax.set_yticklabels([LABELS[k][0] for k in comparisons], fontsize=8.5)
        ax.axvline(x=0, color="#999999", lw=0.8, linestyle="-")
        ax.axvline(x=0.5, color="#999999", lw=0.8, linestyle="--", alpha=0.6)
        ax.text(0.51, len(comparisons) - 0.42, "large effect", fontsize=6.5,
                color="#777777", rotation=90, va="top")
        ax.set_xlabel("Effect size")
        ax.set_title("Effect sizes: top-10 and neural-TF contrasts\n"
                     "circular = score contains the label stream; honest = label streams excluded",
                     fontweight="bold", pad=8)
        ax.legend(fontsize=8, loc="lower right")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        save(fig, "28_effect_sizes")
    except FileNotFoundError as e:
        print(f"  [SKIP] {__file__}: {e}")
        return

if __name__ == "__main__":
    build()
