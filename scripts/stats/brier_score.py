#!/usr/bin/env python
"""Score-separation analysis at the median threshold (Brier context).

2026-09-04 audit fix (Murphy decomposition): the previous "reliability/
resolution" were NOT a valid Murphy decomposition of the Brier score -
the reported resolution (score-scale variance, 0.136) EXCEEDED the
uncertainty (pi*(1-pi), 0.0057) by 24x, which is arithmetically
impossible for a real decomposition (resolution <= uncertainty always).
The score is not a probability, so those keys were mislabeled algebra.

What this version reports:
- the raw Brier score as a DESCRIPTIVE distance between the score scale
  and the binary label (valid arithmetic, no probability claim);
- the classification report at the median threshold (no probability
  assumption);
- a proper Murphy decomposition ONLY as a sanity-check diagnostic on
  binned empirical rates: resolution_b = sum_b (n_b/N)(rate_b - pi)^2,
  reliability_b = sum_b (n_b/N)(rate_b - mean_score_b)^2, with
  Brier = uncertainty - resolution + reliability holding on the binned
  quantities (and resolution <= uncertainty guaranteed by construction).

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

N_BINS = 10


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
    uncertainty = prevalence * (1 - prevalence)

    # ---- Proper Murphy decomposition on binned empirical rates ---------
    # Bins: equal-count deciles of the score. rate_b = positives rate in
    # bin b; mean_score_b = mean score in bin b.
    order = np.argsort(y_score, kind="stable")
    bins = np.array_split(order, N_BINS)
    resolution_b = 0.0   # sum (n_b/N)(rate_b - pi)^2   <= uncertainty
    reliability_b = 0.0  # sum (n_b/N)(rate_b - mean_score_b)^2
    for idx in bins:
        if len(idx) == 0:
            continue
        w_b = len(idx) / n
        rate_b = y_true[idx].mean()
        mean_score_b = y_score[idx].mean()
        resolution_b += w_b * (rate_b - prevalence) ** 2
        reliability_b += w_b * (rate_b - mean_score_b) ** 2
    # Identity check: Brier = uncertainty - resolution + reliability
    # holds exactly only for probability scores; for a non-probability
    # score the identity breaks by the score-scale offset - we report the
    # decomposition of the BINNED rates (valid) and the raw Brier
    # separately, without claiming the identity.
    binned_brier = uncertainty - resolution_b + reliability_b

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
        "brier_score_descriptive": float(brier),
        "brier_interpretation": (
            "mean((score - y)^2) as a scale-distance diagnostic ONLY - the "
            "integrated score is an evidence-weight composite, not a "
            "calibrated probability of RNAi validation"
        ),
        "murphy_decomposition_binned": {
            "n_bins": N_BINS,
            "uncertainty": float(uncertainty),
            "resolution": float(resolution_b),
            "reliability": float(reliability_b),
            "binned_brier_identity": float(binned_brier),
            "note": (
                "Decomposition of binned empirical positive-rates "
                "(deciles of the score). Resolution <= uncertainty holds "
                "by construction. The identity Brier = uncertainty - "
                "resolution + reliability does NOT hold for the raw "
                "score because the score is not a probability."
            ),
        },
        "score_column": score_col,
        "interpretation_note": (
            "The integrated score is an evidence-weight composite, not a "
            "calibrated probability; the Brier number is a descriptive "
            "separation diagnostic. The median-threshold classification "
            "and precision/recall are the valid summaries."
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
    print(f"Brier (descriptive scale-distance): {brier:.6f}")
    print(f"Binned Murphy: uncertainty={uncertainty:.6f} "
          f"resolution={resolution_b:.6f} reliability={reliability_b:.6f} "
          f"(resolution <= uncertainty: {resolution_b <= uncertainty + 1e-12})")
    print(f"Median-threshold classification: acc={accuracy:.3f} "
          f"prec={precision:.3f} rec={recall:.3f} F1={f1:.3f}")

    out_path = RESULTS_DIR / "brier_score.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
