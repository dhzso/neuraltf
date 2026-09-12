"""Precision-Recall (PR) curve with circularity control across benchmarks.

Dual-panel evaluation comparing:
- Panel a: RNAi-screened benchmark (King et al. 2024 mmc5 screening cohort, n=67)
- Panel b: FISH phenotype-confirmed benchmark (King et al. 2024, n=19)

Each panel evaluates:
- Circular benchmark (all 11 streams including RNAi label)
- Honest circularity-controlled model (label-free streams)
- Strict label-free model
- Empirical random background prevalence baseline
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
        precision = np.array(d["pr"]["precision"])
        recall = np.array(d["pr"]["recall"])
        pr_auc = auc(recall, precision)
        ax1.plot(recall, precision, color=color, lw=lw, linestyle=ls,
                 label=f"{label_str} (AP={pr_auc:.3f})")

    base_a = data.get("honest", data.get("circular", {})).get("pr", {}).get("baseline", 0.0057)
    ax1.axhline(y=base_a, color="#888888", lw=0.8, linestyle=":",
                label=f"Baseline ({base_a:.4f})")

    ax1.set_xlabel("Recall (Sensitivity)", fontsize=7.0)
    ax1.set_ylabel("Precision (Positive Predictive Value)", fontsize=7.0)
    ax1.set_title("RNAi-Screened Benchmark (n = 67)", fontsize=7.8, pad=6)
    ax1.legend(loc="upper right", fontsize=5.8, frameon=False)
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
        precision = np.array(d["pr"]["precision"])
        recall = np.array(d["pr"]["recall"])
        pr_auc = auc(recall, precision)
        ax2.plot(recall, precision, color=color, lw=lw, linestyle=ls,
                 label=f"{label_str} (AP={pr_auc:.3f})")

    base_b = pheno_data.get("honest", pheno_data.get("circular", {})).get("pr", {}).get("baseline", 0.0016)
    ax2.axhline(y=base_b, color="#888888", lw=0.8, linestyle=":",
                label=f"Baseline ({base_b:.4f})")

    ax2.set_xlabel("Recall (Sensitivity)", fontsize=7.0)
    ax2.set_ylabel("Precision (Positive Predictive Value)", fontsize=7.0)
    ax2.set_title("FISH Phenotype-Confirmed (n = 19)", fontsize=7.8, pad=6)
    ax2.legend(loc="upper right", fontsize=5.8, frameon=False)
    ax2.set_xlim([-0.02, 1.02])
    ax2.set_ylim([-0.02, 1.02])
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.suptitle("Precision–Recall: Neural TF Recovery Across Benchmarks",
                 fontsize=8.5, fontweight="bold", y=0.99)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.14, wspace=0.22)
    save(fig, "27_pr_curve")


if __name__ == "__main__":
    build()
