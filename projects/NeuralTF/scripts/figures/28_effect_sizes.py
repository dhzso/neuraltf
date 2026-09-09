"""Effect sizes — Cliff's delta + Hedges' g bar panel with honesty labels."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import json

# human-readable labels incl. circularity status
LABELS = {
    "top10_vs_rest": ("Top 10 vs rest", True),
    "neural_vs_non_neural": ("Neural vs rest (all streams)", True),
    "neural_vs_non_neural_honest": ("Neural vs rest (label-free)", False),
    "neural_vs_non_neural_honest_strict": ("Neural vs rest (strict label-free)", False),
}

def build():
    try:
        data_path = RES / "effect_sizes.json"
        with open(data_path) as f:
            data = json.load(f)

        comparisons = [k for k in LABELS if k in data]
        if not comparisons:
            raise ValueError("effect_sizes.json carries no known comparisons")

        fig, ax = plt.subplots(figsize=(W_15COL, 3.4))

        y = np.arange(len(comparisons))
        deltas = [data[k]["cliffs_delta"] for k in comparisons]
        gs = [data[k].get("hedges_g", data[k].get("cohens_d", 0.0)) for k in comparisons]
        pvals = [data[k].get("p_value", None) for k in comparisons]

        width = 0.32
        bars_d = ax.barh(y + width/2, deltas, height=width, color=C_A,
                         edgecolor="none", label="Cliff's δ")
        bars_g = ax.barh(y - width/2, gs, height=width, color=C_B,
                         edgecolor="none", label="Hedges' g")

        for i, (d, g, p) in enumerate(zip(deltas, gs, pvals)):
            ax.text(d + 0.03, y[i] + width/2, f"δ = {d:.2f}", va="center", ha="left",
                    fontsize=6.0, color=C_A)
            if p is not None:
                if p < 1e-15:
                    p_math = r"P < 10^{-15}"
                else:
                    base, exp = f"{p:.1e}".split("e")
                    p_math = rf"P = {base} \times 10^{{{int(exp)}}}"
                g_txt = f"g = {g:.2f} (${p_math}$)"
            else:
                g_txt = f"g = {g:.2f}"
            ax.text(g + 0.03, y[i] - width/2, g_txt, va="center", ha="left",
                    fontsize=6.0, color=C_B)

        ax.set_yticks(y)
        ax.set_yticklabels([LABELS[k][0] for k in comparisons], fontsize=6.8)
        ax.axvline(x=0, color="#555555", lw=0.6)
        ax.axvline(x=0.5, color="#999999", lw=0.6, linestyle=":", label="Reference (0.5)")
        ax.set_xlabel("Effect size (Cliff's δ and Hedges' g)", fontsize=7.0)
        ax.set_title("Effect Sizes and Significance Across Candidate Cohorts",
                     fontsize=8.0, pad=6)
        ax.legend(fontsize=6.2, loc="upper right", frameon=False)
        ax.set_xlim(-0.05, max(max(deltas), max(gs)) * 1.38)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        save(fig, "28_effect_sizes")

    except FileNotFoundError as e:
        print(f"  [SKIP] {__file__}: {e}")
        return

if __name__ == "__main__":
    build()
