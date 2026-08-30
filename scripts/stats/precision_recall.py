#!/usr/bin/env python
"""Precision-recall analysis using RNAi-validated TFs as ground truth.

Computes precision@k, precision@10, and PR-AUC for the integrated
scoring pipeline against known RNAi-validated neural TFs.

Usage:
    python scripts/stats/precision_recall.py
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
FIG_DIR = REPO / "projects" / "NeuralTF" / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


def load_candidates():
    """Load candidate table with scores and proof_status."""
    for fname in ["supplementary_table_S2_fixed_all249.csv", "fstf_ranked_all.csv"]:
        p = RESULTS_DIR / fname
        if p.exists():
            return pd.read_csv(p)
    p = RUN_DIR / "rank.csv"
    if p.exists():
        return pd.read_csv(p)
    raise FileNotFoundError("No candidate score file found")


def compute_pr_curve(y_true, y_scores):
    """Compute precision-recall curve manually."""
    order = np.argsort(-y_scores)
    y_sorted = y_true[order]
    n_pos = y_true.sum()
    precisions = []
    recalls = []
    tp = 0
    for i, label in enumerate(y_sorted):
        if label:
            tp += 1
        precisions.append(tp / (i + 1))
        recalls.append(tp / n_pos if n_pos > 0 else 0)
    return np.array(recalls), np.array(precisions)


def _trapezoid(y, x):
    """Trapezoid rule integration compatible with NumPy <2.0 and >=2.0."""
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    if hasattr(np, "trapz"):
        return float(np.trapz(y, x))
    # Manual fallback
    y = np.asanyarray(y)
    x = np.asanyarray(x)
    return float(np.sum((x[1:] - x[:-1]) * (y[1:] + y[:-1]) / 2.0))


def compute_roc_curve(y_true, y_scores):
    """Compute ROC curve manually without mandatory sklearn dependency."""
    try:
        from sklearn.metrics import roc_curve
        fpr, tpr, _ = roc_curve(y_true, y_scores)
        return fpr, tpr
    except Exception:
        order = np.argsort(-y_scores)
        y_sorted = y_true[order]
        n_pos = int(y_true.sum())
        n_neg = int(len(y_true) - n_pos)
        tpr = [0.0]
        fpr = [0.0]
        tp = fp = 0
        for label in y_sorted:
            if label:
                tp += 1
            else:
                fp += 1
            tpr.append(tp / n_pos if n_pos > 0 else 0.0)
            fpr.append(fp / n_neg if n_neg > 0 else 0.0)
        return np.array(fpr), np.array(tpr)


def main():
    print("=== Precision-Recall Analysis ===")

    df = load_candidates()
    print(f"Loaded {len(df)} candidates")

    gene_col = "gene_id" if "gene_id" in df.columns else "gene_id_v6"
    score_col = None
    for c in ["integrated_score", "composite_score", "dirichlet_median_score", "fixed_weight_score"]:
        if c in df.columns:
            score_col = c
            break
    if score_col is None:
        print("Error: no score column found")
        return 1

    if "proof_status" not in df.columns:
        print("Error: proof_status column not found")
        return 1

    y_true = (df["proof_status"] == "known_rnai_validated").astype(int).values
    y_scores = df[score_col].fillna(0).values

    n_pos = int(y_true.sum())
    n_total = int(len(y_true))
    print(f"Ground truth: {n_pos}/{n_total} RNAi-validated TFs")

    recalls, precisions = compute_pr_curve(y_true, y_scores)
    fpr, tpr = compute_roc_curve(y_true, y_scores)

    ap = float(np.mean(precisions[y_true[np.argsort(-y_scores)] == 1])) if n_pos > 0 else 0.0

    k_values = [5, 10, 15, 20]
    prec_at_k = {}
    order = np.argsort(-y_scores)
    for k in k_values:
        top_k_true = y_true[order[:k]]
        prec_at_k[f"precision@{k}"] = float(top_k_true.sum() / k)
        print(f"  Precision@{k}: {top_k_true.sum()}/{k} = {top_k_true.sum()/k:.4f}")

    if len(recalls) > 1:
        pr_auc = _trapezoid(precisions, recalls)
    else:
        pr_auc = 0.0

    if len(fpr) > 1:
        roc_auc = _trapezoid(tpr, fpr)
    else:
        roc_auc = 0.5

    baseline = float(n_pos / n_total) if n_total > 0 else 0.0

    results = {
        "n_candidates": int(n_total),
        "n_positives": int(n_pos),
        "baseline_rate": float(baseline),
        "average_precision": float(ap),
        "pr_auc": float(pr_auc),
        "roc_auc": float(roc_auc),
        "precision_at_k": prec_at_k,
        "score_column": score_col,
        "roc": {
            "fpr": [float(x) for x in fpr],
            "tpr": [float(x) for x in tpr],
            "auc": float(roc_auc),
        },
        "pr": {
            "recall": [float(x) for x in recalls],
            "precision": [float(x) for x in precisions],
            "auc": float(pr_auc),
            "baseline": float(baseline),
        },
    }
    print(f"  PR-AUC: {pr_auc:.4f}")
    print(f"  ROC-AUC: {roc_auc:.4f}")
    print(f"  Average Precision: {ap:.4f}")
    print(f"  Baseline (prevalence): {baseline:.4f}")

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.plot(recalls, precisions, "b-", linewidth=2, label=f"Integrated scoring (PR-AUC={pr_auc:.3f})")
    ax.axhline(y=baseline, color="gray", linestyle="--", label=f"Baseline = {baseline:.3f}")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Precision-Recall Curve: Neural TF Prediction")
    ax.legend(loc="best")
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1.05])
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig_path = FIG_DIR / "precision_recall_curve.png"
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)
    print(f"Saved: {fig_path}")

    out_path = RESULTS_DIR / "precision_recall.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
