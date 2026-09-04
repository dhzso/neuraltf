#!/usr/bin/env python
"""Score-separation analysis at the median threshold (Brier context).

WS3 fix: the previous version treated the integrated evidence score as a
calibrated probability and reported a "Brier skill score". Because the
score is an evidence-weight composite (not P(RNAi-validated)), the Brier
skill was meaningless (it measured -0.66 purely because scores cluster
above 0 while prevalence is ~0.13 — an artifact of scale, not skill).

What IS valid and is kept:
  - raw Brier score and its decomposition (reliability/resolution/
    uncertainty) as DESCRIPTIVE diagnostics of score separation
  - a median-threshold classification report (accuracy/precision/
    recall/F1) which makes no probability assumption
  - the rank-discrimination interpretation: resolution vs uncertainty

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
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def main():
    print("=== Score-Separation Diagnostics (descriptive Brier context) ===")

    p = RUN_DIR / "rank.csv"
    if not p.exists():
        print("Error: run the pipeline first (rank.csv missing)")
        return 1
    df = pd.read_csv(p).drop_duplicates(subset="gene_id", keep="first")

    score_col = None
    for c in ["integrated_score", "composite_score"]:
        if c in df.columns:
            score_col = c
            break
    if score_col is None:
        print("Error: no score column found")
        return 1

    df = df.dropna(subset=[score_col])
    df["is_positive"] = (df["proof_status"] == "known_rnai_validated").astype(int)

    y_true = df["is_positive"].values.astype(float)
    y_score = df[score_col].values.astype(float)
    n = len(y_true)
    n_pos = y_true.sum()
    prevalence = n_pos / n

    brier = np.mean((y_score - y_true) ** 2)
    brier_baseline = prevalence * (1 - prevalence)
    resolution = np.mean((y_score - prevalence) ** 2)
    reliability = brier - resolution + brier_baseline

    y_pred = (y_score >= np.median(y_score)).astype(float)
    tp = np.sum((y_pred == 1) & (y_true == 1))
    fp = np.sum((y_pred == 1) & (y_true == 0))
    tn = np.sum((y_pred == 0) & (y_true == 0))
    fn = np.sum((y_pred == 0) & (y_true == 1))
    accuracy = (tp + tn) / n if n else 0
    precision = tp / (tp + fp) if (tp + fp) else 0
    recall = tp / (tp + fn) if (tp + fn) else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0

    results = {
        "n_candidates": int(n),
        "n_positives": int(n_pos),
        "prevalence": float(prevalence),
        "brier_score": float(brier),
        "brier_baseline": float(brier_baseline),
        "reliability": float(reliability),
        "resolution": float(resolution),
        "uncertainty": float(brier_baseline),
        "score_column": score_col,
        "interpretation_note": (
            "The integrated score is an evidence-weight composite, not a "
            "calibrated probability; the Brier numbers are descriptive "
            "separation diagnostics. Skill-vs-baseline probability claims "
            "were removed (see WS3 audit). The median-threshold "
            "classification and precision/recall are the valid summaries."
        ),
        "classification_at_median": {
            "accuracy": float(accuracy),
            "precision": float(precision),
            "recall": float(recall),
            "f1": float(f1),
            "tp": int(tp), "fp": int(fp), "tn": int(tn), "fn": int(fn),
        },
    }

    print(f"Candidates: {n}; positives: {int(n_pos)} ({prevalence:.4f})")
    print(f"Brier score (descriptive): {brier:.6f}")
    print(f"Resolution / uncertainty: {resolution:.4f} / {brier_baseline:.4f}")
    print(f"Median-threshold classification: acc={accuracy:.3f} "
          f"prec={precision:.3f} rec={recall:.3f} F1={f1:.3f}")

    out_path = RESULTS_DIR / "brier_score.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
