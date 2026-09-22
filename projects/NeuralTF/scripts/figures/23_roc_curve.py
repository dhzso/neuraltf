"""Receiver Operating Characteristic (ROC) curve with circularity control across benchmarks.

Dual-panel evaluation comparing:
- Panel a: RNAi-screened benchmark (King et al. 2024 mmc5 screening cohort)
- Panel b: FISH phenotype-confirmed benchmark (King et al. 2024)

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

# 2026-09-19: cohort sizes computed from the ground-truth module + rank
# table (previously hard-coded "n = 68"/"n = 19" in titles).
sys.path.insert(0, str(REPO / "src"))
from bioforge.evidence.groundtruth import PHENOTYPE_CONFIRMED_V6  # noqa: E402


def _cohort_sizes():
    all_cand = load_all()
    n_screened = int((all_cand["proof_status"] == "tested").sum())
    n_pheno = int(all_cand["gene_id"].astype(str).isin(PHENOTYPE_CONFIRMED_V6).sum())
    return n_screened, n_pheno

# 2026-09-19: cohort sizes computed from the ground-truth module + rank
# table (previously hard-coded "n = 68"/"n = 19" in titles and docstring).
sys.path.insert(0, str(REPO / "src"))
from bioforge.evidence.groundtruth import PHENOTYPE_CONFIRMED_V6  # noqa: E402


def _cohort_sizes():
    all_cand = load_all()
    n_screened = int((all_cand["proof_status"] == "tested").sum())
    n_pheno = int(all_cand["gene_id"].astype(str).isin(PHENOTYPE_CONFIRMED_V6).sum())
    return n_screened, n_pheno


def build():
    data_path = RES / "precision_recall.json"
    if not data_path.exists():
        raise FileNotFoundError(
            f"{data_path} missing — run scripts/stats/precision_recall.py first"
        )
    with open(data_path) as f:
        data = json.load(f)

    n_screened, n_pheno = _cohort_sizes()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_2COL, 3.4), dpi=500)

    specs = [
        ("circular", "#687787", "All 11 streams (circular)", 1.4, "--"),
        ("honest", C_A, "Circularity-controlled", 1.8, "-"),
        ("honest_strict", C_B, "Strict label-free", 1.5, "-"),
    ]

    # Panel a: RNAi-screened cohort (n = 68)
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
    ax1.set_title(f"RNAi-Screened Benchmark (n = {n_screened})", fontsize=7.8, pad=6)
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
    ax2.set_title(f"FISH Phenotype-Confirmed (n = {n_pheno})", fontsize=7.8, pad=6)
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
