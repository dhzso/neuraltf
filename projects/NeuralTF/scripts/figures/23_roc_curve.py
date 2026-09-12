"""Receiver Operating Characteristic (ROC) curve with circularity control.

Evaluates recovery of RNAi-screened neural regulators (King 2024 mmc5
screening list; phenotype status tracked separately via
`phenotype_confirmed`) comparing:
- Circular benchmark (all 11 streams including RNAi label)
- Honest circularity-controlled model (label-free streams)
- Strict label-free model
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
        fpr = np.array(d["roc"]["fpr"])
        tpr = np.array(d["roc"]["tpr"])
        roc_auc = auc(fpr, tpr)
        ax.plot(fpr, tpr, color=color, lw=lw, linestyle=ls,
                label=f"{label_str} (AUC = {roc_auc:.3f})")

    ax.plot([0, 1], [0, 1], color="#888888", lw=0.8, linestyle=":",
            label="Random classifier (AUC = 0.500)")

    ax.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=7.0)
    ax.set_ylabel("True Positive Rate (Sensitivity)", fontsize=7.0)
    ax.set_title("Receiver Operating Characteristic: Neural TF Recovery",
                 fontsize=8.0, pad=6)
    ax.legend(loc="lower right", fontsize=6.2, frameon=False)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.subplots_adjust(left=0.14, right=0.96, top=0.90, bottom=0.14)
    save(fig, "23_roc_curve")


if __name__ == "__main__":
    build()
