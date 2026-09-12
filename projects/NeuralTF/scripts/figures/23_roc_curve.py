"""Receiver Operating Characteristic (ROC) curve with circularity control across benchmarks.

Dual-panel evaluation comparing:
- Panel a: RNAi-screened benchmark (King et al. 2024 mmc5 screening cohort, n=67)
- Panel b: FISH phenotype-confirmed benchmark (King et al. 2024, n=19)

Each panel evaluates:
- Circular benchmark (all 11 streams including RNAi label)
- Honest circularity-controlled model (label-free streams)
- Strict label-free model
- Random classifier baseline (AUC = 0.500)
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

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_2COL, 3.4), dpi=500)

    specs = [
        ("circular", "#687787", "All 11 streams (circular)", 1.4, "--"),
        ("honest", C_A, "Circularity-controlled", 1.8, "-"),
        ("honest_strict", C_B, "Strict label-free", 1.5, "-"),
    ]

    # Panel a: RNAi-screened cohort (n = 67)
    panel_tag(ax1, "a", x=-0.14, y=1.04)
    for key, color, label_str, lw, ls in specs:
        if key not in data:
            continue
        d = data[key]
        fpr = np.array(d["roc"]["fpr"])
        tpr = np.array(d["roc"]["tpr"])
        roc_auc = auc(fpr, tpr)
        ax1.plot(fpr, tpr, color=color, lw=lw, linestyle=ls,
                 label=f"{label_str} ({roc_auc:.3f})")

    ax1.plot([0, 1], [0, 1], color="#888888", lw=0.8, linestyle=":",
             label="Random (0.500)")
    ax1.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=7.0)
    ax1.set_ylabel("True Positive Rate (Sensitivity)", fontsize=7.0)
    ax1.set_title("RNAi-Screened Benchmark (n = 67)", fontsize=7.8, pad=6)
    ax1.legend(loc="lower right", fontsize=5.8, frameon=False)
    ax1.set_xlim([-0.02, 1.02])
    ax1.set_ylim([-0.02, 1.02])
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel b: FISH phenotype-confirmed cohort (n = 19)
    panel_tag(ax2, "b", x=-0.14, y=1.04)
    pheno_data = data.get("phenotype_confirmed", {})
    for key, color, label_str, lw, ls in specs:
        if key not in pheno_data:
            continue
        d = pheno_data[key]
        fpr = np.array(d["roc"]["fpr"])
        tpr = np.array(d["roc"]["tpr"])
        roc_auc = auc(fpr, tpr)
        ax2.plot(fpr, tpr, color=color, lw=lw, linestyle=ls,
                 label=f"{label_str} ({roc_auc:.3f})")

    ax2.plot([0, 1], [0, 1], color="#888888", lw=0.8, linestyle=":",
             label="Random (0.500)")
    ax2.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=7.0)
    ax2.set_ylabel("True Positive Rate (Sensitivity)", fontsize=7.0)
    ax2.set_title("FISH Phenotype-Confirmed (n = 19)", fontsize=7.8, pad=6)
    ax2.legend(loc="lower right", fontsize=5.8, frameon=False)
    ax2.set_xlim([-0.02, 1.02])
    ax2.set_ylim([-0.02, 1.02])
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.suptitle("Receiver Operating Characteristic: Neural TF Recovery Across Benchmarks",
                 fontsize=8.5, fontweight="bold", y=0.99)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.14, wspace=0.22)
    save(fig, "23_roc_curve")


if __name__ == "__main__":
    build()
