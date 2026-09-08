"""ROC and PR curves using RNAi TFs as ground truth (circularity-controlled).

Reads the dual-evaluation JSON written by scripts/stats/precision_recall.py:
  - "circular": all 9 streams (rnai/neural_* included — inflated)
  - "honest" : label-encoding streams excluded (publishable estimate)
Both ROC and PR panels overlay the two curves.
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
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

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_2COL, 3.4))

    specs = [("circular", "#455A64", "All 9 streams (circular benchmark)"),
             ("honest", C_A, "Label-free (honest lower bound)"),
             ("honest_strict", C_B, "Strict (excl. reproducibility)")]
    for key, color, label in specs:
        if key not in data:
            continue
        d = data[key]
        fpr = np.array(d["roc"]["fpr"])
        tpr = np.array(d["roc"]["tpr"])
        roc_auc = auc(fpr, tpr)
        ax1.plot(fpr, tpr, color=color, lw=1.6,
                 label=f"{label} (AUC = {roc_auc:.3f})")

        precision = np.array(d["pr"]["precision"])
        recall = np.array(d["pr"]["recall"])
        pr_auc = auc(recall, precision)
        ax2.plot(recall, precision, color=color, lw=1.6,
                 label=f"{label} (AP = {pr_auc:.3f})")

    ax1.plot([0, 1], [0, 1], color="#999999", lw=0.8, linestyle="--",
             label="Chance diagonal (AUC = 0.500)")
    ax1.set_xlabel("False positive rate", fontsize=8)
    ax1.set_ylabel("True positive rate (sensitivity)", fontsize=8)
    ax1.set_title("ROC: recovery of RNAi-validated neural TFs",
                  fontweight="bold", fontsize=8.5, pad=6)
    ax1.legend(loc="lower right", fontsize=6.8, frameon=False)
    ax1.set_xlim([-0.02, 1.02])
    ax1.set_ylim([-0.02, 1.02])
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    panel_tag(ax1, "a")

    baseline = data.get("honest", data.get("circular", {})).get(
        "pr", {}).get("baseline", 0.15)
    ax2.axhline(y=baseline, color="#999999", lw=0.8, linestyle="--",
                label=f"Cohort prevalence ({baseline:.2f})")
    ax2.set_xlabel("Recall (sensitivity)", fontsize=8)
    ax2.set_ylabel("Precision (positive predictive value)", fontsize=8)
    ax2.set_title("Precision–Recall: recovery of RNAi-validated TFs",
                  fontweight="bold", fontsize=8.5, pad=6)
    ax2.legend(loc="upper right", fontsize=6.8, frameon=False)
    ax2.set_xlim([-0.02, 1.02])
    ax2.set_ylim([-0.02, 1.02])
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    panel_tag(ax2, "b")

    fig.tight_layout()
    save(fig, "23_roc_pr_curve")


if __name__ == "__main__":
    build()
