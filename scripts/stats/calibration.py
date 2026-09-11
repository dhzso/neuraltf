#!/usr/bin/env python
"""Calibration / discrimination analysis for integrated scores.

Bins integrated scores into deciles and computes the empirical positive
rate (RNAi-validated TFs) per bin. NOTE (WS3): the integrated score is
an evidence-weight score in [0,1], not a calibrated probability, so a
classic "perfect calibration" diagonal is conceptually invalid. We
report two honest discrimination metrics:

  - rank-discrimination error: mean |empirical positive rate - prevalence|
    per decile (how far decile rates deviate from the base rate — this is
    what the old code called "ECE")
  - true ECE against prevalence-weighted decile means (documented as a
    discrimination proxy)

The reliability plot itself (score decile vs observed rate) is valid and
is what figure 30 renders.

Usage:
    python scripts/stats/calibration.py --n-bins 10
"""

import argparse
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


def load_all_genes():
    """Load candidate table."""
    p = RUN_DIR / "rank.csv"
    if p.exists():
        return pd.read_csv(p).drop_duplicates(subset="gene_id", keep="first")
    raise FileNotFoundError("No candidate score file found in runs/pipeline_run/rank.csv")


def main():
    parser = argparse.ArgumentParser(description="Calibration analysis for integrated scores")
    parser.add_argument("--n-bins", type=int, default=10, help="Number of score bins (default 10)")
    args = parser.parse_args()

    print("=== Calibration Analysis ===")

    df = load_all_genes()
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

    df = df.dropna(subset=[score_col]).copy()
    df["is_positive"] = (df["proof_status"] == "tested").astype(int)

    n_pos = df["is_positive"].sum()
    n_total = len(df)
    prevalence = n_pos / n_total if n_total > 0 else 0
    print(f"Candidates: {n_total}, Positives: {n_pos}, Prevalence: {prevalence:.4f}")

    # 2026-09-06 audit note: qcut over rank(method='first') splits tied
    # scores across decile boundaries by row order (verified: the
    # 541-gene tie group at 0.5667 straddles two deciles). With
    # massive score ties this is unavoidable for equal-count binning;
    # the per-decile rates near tie boundaries must be read as
    # approximate. Documented rather than "fixed" (any alternative —
    # score-value binning — makes bins unequal and breaks the
    # equal-count enrichment logic above).
    df["decile"] = pd.qcut(df[score_col].rank(method="first"), q=args.n_bins, labels=False)
    df["decile"] = args.n_bins - 1 - df["decile"]

    bin_stats = []
    for decile in sorted(df["decile"].unique()):
        bin_df = df[df["decile"] == decile]
        n_bin = len(bin_df)
        n_pos_bin = bin_df["is_positive"].sum()
        empirical_rate = n_pos_bin / n_bin if n_bin > 0 else 0
        mean_score = bin_df[score_col].mean()
        lo = bin_df[score_col].min()
        hi = bin_df[score_col].max()

        bin_stats.append({
            "decile": int(decile),
            "score_range_lo": float(lo),
            "score_range_hi": float(hi),
            "mean_score": float(mean_score),
            "n_candidates": int(n_bin),
            "n_positives": int(n_pos_bin),
            "empirical_positive_rate": float(empirical_rate),
        })
        print(f"  Decile {decile}: [{lo:.3f}, {hi:.3f}], n={n_bin}, "
              f"positives={n_pos_bin}, rate={empirical_rate:.4f}")

    bin_df_out = pd.DataFrame(bin_stats)

    mean_scores = bin_df_out["mean_score"].values
    emp_rates = bin_df_out["empirical_positive_rate"].values
    # 2026-09-06 audit fix: the old "calibration error" metrics were
    # construction artifacts. With prevalence pi ~ 0.0057, the TOP decile
    # rate is capped at ~10*pi and the other nine at >= 0; a perfectly
    # ranking score therefore ALWAYS produces max deviation ~= 10*pi -
    # pi = the ceiling, i.e. the metric is a monotone transform of rank
    # discrimination (already captured by AUC) and CANNOT be small for
    # any top-concentrating score. The reported metric is now the
    # TOP-DECILE ENRICHMENT (the quantity the funnel actually needs):
    # fold-enrichment of positives in the top decile over prevalence,
    # with an exact binomial CI and a one-sided p-value. Per-decile
    # binomial CIs are exported for the figure.
    from scipy import stats as _st

    # 2026-09-06 fix: the TOP decile is decile 0 (the code re-labels so
    # decile 0 = highest scores: `args.n_bins - 1 - qcut_labels`). The
    # first version of this fix read bin_stats[-1] — the LOWEST decile —
    # producing "enrichment 0.0x, p=1.0" (64/67 positives live in
    # decile 0, as the bin table shows).
    top = bin_stats[0]
    k_top, n_top = top["n_positives"], top["n_candidates"]
    # exact Clopper-Pearson interval on the top-decile positive rate
    ci_lo = float(_st.beta.ppf(0.025, k_top, n_top - k_top + 1)) if k_top > 0 else 0.0
    ci_hi = float(_st.beta.ppf(0.975, k_top + 1, n_top - k_top)) if k_top < n_top else 1.0
    top_enrichment = (k_top / n_top) / prevalence if prevalence > 0 else float("nan")
    top_p = float(_st.binomtest(k_top, n_top, prevalence,
                                alternative="greater").pvalue) \
        if hasattr(_st, "binomtest") else float(_st.binom_test(k_top, n_top, prevalence))

    per_decile_ci = []
    for b in bin_stats:
        k, nb = b["n_positives"], b["n_candidates"]
        lo_b = float(_st.beta.ppf(0.025, k, nb - k + 1)) if k > 0 else 0.0
        hi_b = float(_st.beta.ppf(0.975, k + 1, nb - k)) if k < nb else 1.0
        per_decile_ci.append({"decile": b["decile"], "rate_lo": lo_b, "rate_hi": hi_b})

    results = {
        "n_bins": args.n_bins,
        "n_total": int(n_total),
        "n_positives": int(n_pos),
        "prevalence": float(prevalence),
        "top_decile_enrichment": {
            "fold_enrichment": float(top_enrichment),
            "top_positive_rate": float(k_top / n_top) if n_top else 0.0,
            "rate_ci95": [ci_lo, ci_hi],
            "p_one_sided_binomial": top_p,
            "n_positives_top_decile": int(k_top),
            "n_candidates_top_decile": int(n_top),
        },
        "retired_metric_note": (
            "mean/max 'calibration error' vs prevalence RETIRED (2026-09-06): "
            "with prevalence pi, a perfectly-ranking score forces max "
            "deviation to ~9*pi by construction (top decile caps at "
            "10*pi, the rest sit at 0) — the metric was a monotone "
            "transform of rank discrimination, not calibration, and "
            "could never be small for any top-concentrating score."
        ),
        "positive_label_note": (
            "positives = proof_status == tested ONLY; the "
            "near-positive known_fstf group (n~61) counts as "
            "negative here, so the metric measures RNAi-validated "
            "discrimination, not general 'neural TF-ness'."
        ),
        "discrimination_note": (
            "integrated_score is an evidence-weight score, not a probability; "
            "the top-decile enrichment (with exact binomial inference) is the "
            "valid discrimination summary, not a reliability diagonal"
        ),
        "bin_centers": [float(b["mean_score"]) for b in bin_stats],
        "observed_fractions": [float(b["empirical_positive_rate"]) for b in bin_stats],
        "per_decile_ci95": per_decile_ci,
        "bin_mean_scores": [float(b["mean_score"]) for b in bin_stats],
        "bin_mean_scores_note": (
            "mean score per decile, for context only - the score is NOT a "
            "probability, so these must not be plotted as an expected/"
            "reliability diagonal against observed fractions"
        ),
        "bin_counts": [int(b["n_candidates"]) for b in bin_stats],
        "bin_stats": bin_stats,
    }
    print(f"\nTop-decile enrichment: {top_enrichment:.2f}x "
          f"(rate {k_top}/{n_top} = {k_top / n_top:.4f} vs prevalence "
          f"{prevalence:.4f}), exact one-sided binomial p = {top_p:.3e}")
    print("Retired mean/max 'calibration error' (construction-bounded "
          "metric — see retired_metric_note).")

    out_path = RESULTS_DIR / "calibration_stats.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
