"""Effect sizes — box/violin with Cliff's delta annotations."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import json

def build():
    data_path = RES / "effect_sizes.json"
    with open(data_path) as f:
        data = json.load(f)

    fig, ax = plt.subplots(figsize=(7, 5))

    comparisons = list(data.keys())
    deltas = [data[k]["cliffs_delta"] for k in comparisons]
    pvals = [data[k].get("p_value", None) for k in comparisons]
    sizes = [abs(d) for d in deltas]
    colors = [C_HL if s > 0.5 else C_A if s > 0.3 else C_B for s in sizes]

    y = np.arange(len(comparisons))
    bars = ax.barh(y, deltas, color=colors, edgecolor="white", lw=0.5, height=0.6)

    for i, (d, p) in enumerate(zip(deltas, pvals)):
        label = f"δ = {d:.3f}"
        if p is not None:
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            label += f" ({sig})"
        ax.text(d + (0.02 if d >= 0 else -0.02), i, label,
                va="center", ha="left" if d >= 0 else "right", fontsize=7, fontweight="bold")

    ax.set_yticks(y)
    ax.set_yticklabels([k.replace("_", " ").title() for k in comparisons], fontsize=8)
    ax.axvline(x=0, color="#999999", lw=0.8, linestyle="-")
    ax.set_xlabel("Cliff's delta (effect size)")
    ax.set_title("Effect sizes between neural TFs and control groups\n"
                 "Cliff's δ > 0.3 = medium effect; > 0.5 = large effect",
                 fontweight="bold", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "28_effect_sizes")

if __name__ == "__main__":
    build()
