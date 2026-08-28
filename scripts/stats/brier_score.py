#!/usr/bin/env python
"""Brier score for integrated scoring as a probabilistic classifier.

Computes the Brier score (mean squared error of probability predictions)
for the integrated scoring pipeline treating RNAi-validated TFs as the
binary ground truth.

Usage:
    python scripts/stats/brier_score.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
FIG_DIR = REPO / "projects" / "NeuralTF" / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


def main():
    print("=== Brier Score Analysis ===")

    for fname in ["supplementary_table_S2_fixed_all249.csv", "fstf_ranked_all.csv"]:
        p = RESULTS_DIR / fname
        if p.exists():
            df = pd.read_csv(p)
            break
    else:
        p = RUN_DIR / "rank.csv"
        if p.exists():
            df = pd.read_csv(p)
        else:
            print("Error: no candidate score file found")
            return 1

    score_col = None
    for c in ["integrated_score", "composite_score", "dirichlet_median_score", "fixed_weight_score"]:
        if c in df.columns:
            score_col = c
            break
    if score_col is None:
        print("Error: no score column found")
        return 1

    gene_col = "gene_id" if "gene_id" in df.columns else "gene_id_v6"
    df = df.dropna(subset=[score_col])
    df["is_positive"] = (df["proof_status"] == "known_rnai_validated").astype(int)

    y_true = df["is_positive"].values.astype(float)
    y_prob = df[score_col].values.astype(float)
    n = len(y_true)
    n_pos = y_true.sum()
    prevalence = n_pos / n

    brier = np.mean((y_prob - y_true) ** 2)
    brier_normalized = brier / (prevalence * (1 - prevalence)) if prevalence > 0 and prevalence < 1 else 0.0

    brier_baseline = prevalence * (1 - prevalence)
    brier_skill = 1.0 - (brier / brier_baseline) if brier_baseline > 0 else 0.0

    resolution = np.mean((y_prob - prevalence) ** 2)
    uncertainty = prevalence * (1 - prevalence)
    reliability = brier - resolution + uncertainty

    y_pred_binary = (y_prob >= np.median(y_prob)).astype(float)
    tp = np.sum((y_pred_binary == 1) & (y_true == 1))
    fp = np.sum((y_pred_binary == 1) & (y_true == 0))
    tn = np.sum((y_pred_binary == 0) & (y_true == 0))
    fn = np.sum((y_pred_binary == 0) & (y_true == 1))

    accuracy = (tp + tn) / n if n > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    results = {
        "n_candidates": int(n),
        "n_positives": int(n_pos),
        "prevalence": float(prevalence),
        "brier_score": float(brier),
        "brier_baseline": float(brier_baseline),
        "brier_skill_score": float(brier_skill),
        "brier_normalized": float(brier_normalized),
        "reliability": float(reliability),
        "resolution": float(resolution),
        "uncertainty": float(uncertainty),
        "score_column": score_col,
        "classification_at_median": {
            "accuracy": float(accuracy),
            "precision": float(precision),
            "recall": float(recall),
            "f1": float(f1),
            "tp": int(tp),
            "fp": int(fp),
            "tn": int(tn),
            "fn": int(fn),
        },
    }

    print(f"Candidates: {n}")
    print(f"Positives: {n_pos} ({prevalence:.4f})")
    print(f"\nBrier score:          {brier:.6f}")
    print(f"Brier baseline:       {brier_baseline:.6f}")
    print(f"Brier skill score:    {brier_skill:.6f}")
    print(f"Reliability (cal.):   {reliability:.6f}")
    print(f"Resolution:           {resolution:.6f}")
    print(f"Uncertainty:          {uncertainty:.6f}")
    print(f"\nClassification at median threshold:")
    print(f"  Accuracy:  {accuracy:.4f}")
    print(f"  Precision: {precision:.4f}")
    print(f"  Recall:    {recall:.4f}")
    print(f"  F1:        {f1:.4f}")

    print(f"\nInterpretation:")
    if brier_skill > 0.5:
        print(f"  Good skill (skill={brier_skill:.3f} > 0.5)")
    elif brier_skill > 0:
        print(f"  Moderate skill (0 < skill={brier_skill:.3f} < 0.5)")
    else:
        print(f"  Poor skill (skill={brier_skill:.3f} <= 0)")

    out_path = RESULTS_DIR / "brier_score.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
