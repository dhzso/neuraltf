"""Precision-Recall (PR) curve with circularity control across benchmarks.

Dual-panel evaluation comparing:
- Panel a: RNAi-screened benchmark (King et al. 2024 mmc5 screening cohort)
- Panel b: FISH phenotype-confirmed benchmark (King et al. 2024)

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

# 2026-09-19: cohort sizes + PR baselines computed from the ground-truth
# module + rank table (previously hard-coded titles n=68/n=19 and fallback
# baselines 0.0057/0.0016 that understated the true prevalence).
sys.path.insert(0, str(REPO / "src"))
from bioforge.evidence.groundtruth import PHENOTYPE_CONFIRMED_V6  # noqa: E402


def _cohort_sizes():
    all_cand = load_all()
    n_screened = int((all_cand["proof_status"] == "tested").sum())
    n_pheno = int(all_cand["gene_id"].astype(str).isin(PHENOTYPE_CONFIRMED_V6).sum())
    n_total = int(len(all_cand))
    return n_screened, n_pheno, n_total


def build():
    data_path = RES / "precision_recall.json"
    if not data_path.exists():
        raise FileNotFoundError(
            f"{data_path} missing — run scripts/stats/precision_recall.py first"
        )
    with open(data_path) as f:
        data = json.load(f)

    n_screened, n_pheno, n_total = _cohort_sizes()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_2COL, 3.4), dpi=500)

    specs = [
        ("circular", "#687787", "All 11 streams (circular)", 1.4, "--"),
        ("honest", C_A, "Circularity-controlled (label-free)", 1.8, "-"),
        ("honest_strict", C_B, "Strict label-free", 1.5, "-"),
    ]

    # Panel a: RNAi-screened cohort (n = 68)
    panel_tag(ax1, "a", x=-0.14, y=1.04)
    for key, color, label_str, lw, ls in specs:
        if key not in data:
            continue
        d = data[key]
        precision = np.array(d["pr"]["precision"])
        recall = np.array(d["pr"]["recall"])
        # 2026-09-13: trapezoidal area over the PR curve — labelled as
        # PR-AUC (NOT "AP"; average precision is a different, non-
        # interpolated quantity and the two differ, e.g. 0.047 vs 0.049).
        pr_auc = auc(recall, precision)
        ax1.plot(recall, precision, color=color, lw=lw, linestyle=ls,
                 label=f"{label_str} (PR-AUC={pr_auc:.3f})")

    base_a = data.get("honest", data.get("circular", {})).get("pr", {}).get(
        "baseline", n_screened / n_total)
    ax1.axhline(y=base_a, color="#888888", lw=0.8, linestyle=":",
                label=f"Prevalence baseline ({base_a:.4f})")

    ax1.set_xlabel("Recall (Sensitivity)", fontsize=7.0)
    ax1.set_ylabel("Precision (Positive Predictive Value)", fontsize=7.0)
    ax1.set_title(f"RNAi-Screened Benchmark (n = {n_screened})", fontsize=7.8, pad=6)
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
                 label=f"{label_str} (PR-AUC={pr_auc:.3f})")

    base_b = pheno_data.get("honest", pheno_data.get("circular", {})).get(
            "pr", {}).get("baseline", n_pheno / n_total)
    ax2.axhline(y=base_b, color="#888888", lw=0.8, linestyle=":",
                label=f"Prevalence baseline ({base_b:.4f})")

    ax2.set_xlabel("Recall (Sensitivity)", fontsize=7.0)
    ax2.set_ylabel("Precision (Positive Predictive Value)", fontsize=7.0)
    ax2.set_title(f"FISH Phenotype-Confirmed (n = {n_pheno})", fontsize=7.8, pad=6)
    ax2.legend(loc="upper right", fontsize=5.8, frameon=False)
    ax2.set_xlim([-0.02, 1.02])
    ax2.set_ylim([-0.02, 1.02])
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # 2026-09-23: header raised ~2.5 pt + subtitle pulled 3.5 pt closer to
    # the title so the subtitle clears the panel tag "b" and panel titles
    # below it (previously the bboxes overlapped by ~1.7 pt).
    title_block(
        fig,
        "Precision–Recall: Neural TF Recovery Across Benchmarks",
        f"Benchmarks: RNAi-screened ($n$ = {n_screened}) and FISH-confirmed ($n$ = {n_pheno}) cohorts",
        y=0.9952, sub_y=0.9434,
    )
    fig.subplots_adjust(left=0.08, right=0.98, top=0.84, bottom=0.16, wspace=0.22)
    save(fig, "27_pr_curve")


if __name__ == "__main__":
    build()
