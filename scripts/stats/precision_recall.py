#!/usr/bin/env python
"""Precision-recall analysis using RNAi-validated TFs as ground truth.

WS3 fix for circularity: the `rnai` and `neural_enriched` streams ARE
part of the integrated score, and the ground truth (proof_status ==
known_rnai_validated) is derived from the same King mmc5 RNAi table that
feeds the rnai stream — so a naive evaluation is circular and inflates
ROC-AUC (~0.91 vs ~0.69 honest). This script therefore reports BOTH:

  1. circular   — score includes rnai + neural_enriched streams
  2. honest     — score recomputed EXCLUDING rnai, neural_enriched, and
                  neural_specificity (any stream directly encoding the
                  RNAi/neural labels)

The honest number is the publishable estimate of how well the remaining
multi-atlas evidence recovers known neural TFs.

Outputs:
  results/precision_recall.json (curves for both evaluations)
  figure rendered by the numbered script (23_roc_pr_curve.py).

Usage:
    python scripts/stats/precision_recall.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

STREAMS = ["expression", "specificity", "reproducibility", "rnai",
           "correlation", "neural_enriched", "neural_specificity",
           "perez_lineage", "perez_influence"]
W_DEFAULT = {"expression": 0.2, "specificity": 0.1, "reproducibility": 0.1,
             "rnai": 0.1, "correlation": 0.1, "neural_enriched": 0.1,
             "neural_specificity": 0.1, "perez_lineage": 0.1,
             "perez_influence": 0.1}
# streams that directly encode the ground-truth labels
CIRCULAR_STREAMS = {"rnai", "neural_enriched", "neural_specificity"}


def load_candidates():
    p = RUN_DIR / "rank.csv"
    if not p.exists():
        raise FileNotFoundError(f"No candidate score file at {p}")
    return pd.read_csv(p).drop_duplicates(subset="gene_id", keep="first")


def recompute_excluding_circular(df: pd.DataFrame) -> np.ndarray:
    """Renormalized weighted score over the non-circular streams only."""
    keep = [s for s in STREAMS if s in df.columns and s not in CIRCULAR_STREAMS]
    S = df[keep].to_numpy(dtype=float)
    W = np.array([W_DEFAULT[s] for s in keep])
    valid = ~np.isnan(S)
    S_filled = np.nan_to_num(S, nan=0.0)
    num = S_filled @ W
    den = valid.astype(float) @ W
    den = np.where(den > 0, den, 1.0)
    return num / den


def compute_pr_curve(y_true, y_scores):
    order = np.argsort(-y_scores, kind="stable")
    y_sorted = y_true[order]
    n_pos = y_true.sum()
    precisions, recalls = [], []
    tp = 0
    for i, label in enumerate(y_sorted):
        if label:
            tp += 1
        precisions.append(tp / (i + 1))
        recalls.append(tp / n_pos if n_pos > 0 else 0)
    return np.array(recalls), np.array(precisions)


def _trapezoid(y, x):
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    if hasattr(np, "trapz"):
        return float(np.trapz(y, x))
    y = np.asanyarray(y); x = np.asanyarray(x)
    return float(np.sum((x[1:] - x[:-1]) * (y[1:] + y[:-1]) / 2.0))


def compute_roc_curve(y_true, y_scores):
    order = np.argsort(-y_scores, kind="stable")
    y_sorted = y_true[order]
    n_pos = int(y_true.sum())
    n_neg = int(len(y_true) - n_pos)
    tpr, fpr = [0.0], [0.0]
    tp = fp = 0
    for label in y_sorted:
        if label:
            tp += 1
        else:
            fp += 1
        tpr.append(tp / n_pos if n_pos > 0 else 0.0)
        fpr.append(fp / n_neg if n_neg > 0 else 0.0)
    return np.array(fpr), np.array(tpr)


def evaluate(df, score_col, label):
    y_true = (df["proof_status"] == "known_rnai_validated").astype(int).values
    y_scores = pd.to_numeric(df[score_col], errors="coerce").fillna(0).values
    n_pos = int(y_true.sum())

    recalls, precisions = compute_pr_curve(y_true, y_scores)
    fpr, tpr = compute_roc_curve(y_true, y_scores)
    ap = float(np.mean(precisions[y_true[np.argsort(-y_scores, kind="stable")] == 1])) \
        if n_pos > 0 else 0.0
    pr_auc = _trapezoid(precisions, recalls) if len(recalls) > 1 else 0.0
    roc_auc = _trapezoid(tpr, fpr) if len(fpr) > 1 else 0.5

    prec_at_k = {}
    order = np.argsort(-y_scores, kind="stable")
    for k in (5, 10, 15, 20):
        top_k_true = y_true[order[:k]]
        prec_at_k[f"precision@{k}"] = float(top_k_true.sum() / k)

    print(f"\n[{label}] n={len(y_true)}, positives={n_pos}")
    print(f"  PR-AUC: {pr_auc:.4f}  ROC-AUC: {roc_auc:.4f}  AP: {ap:.4f}")
    for k, v in prec_at_k.items():
        print(f"  {k}: {v:.4f}")

    return {
        "n_candidates": int(len(y_true)),
        "n_positives": n_pos,
        "average_precision": ap,
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
            "baseline": float(n_pos / len(y_true)) if len(y_true) else 0.0,
        },
    }


def main():
    print("=== Precision-Recall Analysis (circularity-controlled) ===")

    df = load_candidates()
    print(f"Loaded {len(df)} unique candidates")

    if "proof_status" not in df.columns:
        print("Error: proof_status column not found")
        return 1

    # 1) Circular evaluation (as previously reported)
    circular = evaluate(df, "integrated_score", "circular (all 9 streams)")

    # 2) Honest evaluation (label-encoding streams excluded)
    df["_honest_score"] = recompute_excluding_circular(df)
    honest = evaluate(df, "_honest_score",
                     "honest (rnai/neural_* streams excluded)")

    results = {
        "circular": circular,
        "honest": honest,
        "circularity_note": (
            "The 'circular' evaluation includes the rnai/neural_enriched/"
            "neural_specificity streams that share the King mmc5 ground truth; "
            "the 'honest' evaluation recomputes the renormalized weighted "
            "score excluding those streams and is the publishable estimate."
        ),
    }

    out_path = RESULTS_DIR / "precision_recall.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
