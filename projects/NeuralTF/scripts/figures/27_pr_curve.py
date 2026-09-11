"""Precision-Recall (PR) curve with circularity control.

Evaluates precision and average precision (AP) across recall thresholds for neural TF recovery:
- Circular benchmark (all 9 streams)
- Honest circularity-controlled model (label-free streams)
- Strict label-free model
- Empirical background prevalence baseline
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import json
from sklearn.metrics import auc


def build():
    data_path = RES / "precision_recall.json"
    if not data_path.exists():
        raise FileNotFoundError(
            f"{data_path} missing — run scripts/stats/precision_recall.py first"
        )
    with open(data_path) as f:
        data = json.load(f)

    fig, ax = plt.subplots(figsize=(W_15COL, 3.6), dpi=500)

    specs = [
        ("circular", "#687787", "All 11 streams (circular)", 1.4, "--"),
        ("honest", C_A, "Circularity-controlled (label-free)", 1.8, "-"),
        ("honest_strict", C_B, "Strict label-free", 1.5, "-"),
    ]

    for key, color, label_str, lw, ls in specs:
        if key not in data:
            continue
        d = data[key]
        precision = np.array(d["pr"]["precision"])
        recall = np.array(d["pr"]["recall"])
        pr_auc = auc(recall, precision)
        ax.plot(recall, precision, color=color, lw=lw, linestyle=ls,
                label=f"{label_str} (AP = {pr_auc:.3f})")

    baseline = data.get("honest", data.get("circular", {})).get(
        "pr", {}).get("baseline", 0.15)
    ax.axhline(y=baseline, color="#888888", lw=0.8, linestyle=":",
               label=f"Random baseline ({baseline:.2f})")

    ax.set_xlabel("Recall (Sensitivity)", fontsize=7.0)
    ax.set_ylabel("Precision (Positive Predictive Value)", fontsize=7.0)
    ax.set_title("Precision–Recall: Neural TF Recovery",
                 fontsize=8.0, pad=26)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=2, frameon=False,
              fontsize=5.8, handletextpad=0.3, columnspacing=1.0)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.subplots_adjust(left=0.14, right=0.96, top=0.82, bottom=0.14)
    save(fig, "27_pr_curve")


if __name__ == "__main__":
    build()
