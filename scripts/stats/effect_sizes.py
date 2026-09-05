#!/usr/bin/env python
"""Effect size analysis for score distributions.

Computes Cliff's delta, Cohen's d, and Hedges' g for the key contrasts,
plus Mann-Whitney U tests.

2026-09-04 audit notes (documented circularity + corrections):
- Top-10 vs rest is TAUTOLOGICAL (groups defined by the tested score) -
  reported for continuity, flagged as such.
- Neural (RNAi-validated) vs rest on the raw integrated score is
  CIRCULAR (the rnai stream IS the label): the honest variant recomputes
  the score excluding rnai/neural_enriched/neural_specificity.
- Hedges' g (small-sample bias-corrected d, J = 1 - 3/(4(n1+n2)-9)) is
  added alongside Cohen's d; at n1=10 the correction is ~4%.

Usage:
    python scripts/stats/effect_sizes.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
FIG_DIR = REPO / "projects" / "NeuralTF" / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

STREAMS = ["expression", "specificity", "reproducibility", "rnai",
           "correlation", "neural_enriched", "neural_specificity",
           "perez_lineage", "perez_influence"]
W_DEFAULT = np.array([0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
HONEST_EXCLUDE = ["rnai", "neural_enriched", "neural_specificity"]


def cliffs_delta(x, y):
    """Compute Cliff's delta effect size. Returns delta in [-1, 1]."""
    n_x, n_y = len(x), len(y)
    if n_x == 0 or n_y == 0:
        return 0.0
    dominance = 0
    for xi in x:
        dominance += np.sum(xi > y) - np.sum(xi < y)
    return dominance / (n_x * n_y)


def cohens_d(x, y):
    """Compute Cohen's d pooled effect size (ddof=1 sample SDs)."""
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return 0.0
    pooled_std = np.sqrt(((nx - 1) * np.var(x, ddof=1) + (ny - 1) * np.var(y, ddof=1)) / (nx + ny - 2))
    if pooled_std == 0:
        return 0.0
    return (np.mean(x) - np.mean(y)) / pooled_std


def hedges_g(x, y):
    """Hedges' g: Cohen's d times the small-sample bias correction
    J = 1 - 3/(4*(n1+n2) - 9)."""
    d = cohens_d(x, y)
    n = len(x) + len(y)
    if n < 3:
        return d
    return d * (1.0 - 3.0 / (4.0 * n - 9.0))


def honest_score(df: pd.DataFrame) -> pd.Series:
    """Renormalized weighted score over label-independent streams only
    (identical to EvidenceScorer's per-record renormalization)."""
    keep = [s for s in STREAMS if s in df.columns and s not in HONEST_EXCLUDE]
    W = np.array([W_DEFAULT[STREAMS.index(s)] for s in keep])
    M = df[keep].to_numpy(dtype=float)
    mask = ~np.isnan(M)
    num = np.where(mask, M, 0.0) @ W
    den = (mask * W).sum(axis=1)
    return pd.Series(np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0),
                     index=df.index)


def main():
    print("=== Effect Size Analysis ===")

    p = RUN_DIR / "rank.csv"
    if not p.exists():
        print("Error: run the pipeline first (rank.csv missing)")
        return 1
    df = pd.read_csv(p).drop_duplicates(subset="gene_id", keep="first")

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

    df_sorted = df.sort_values(score_col, ascending=False)
    top10_genes = set(df_sorted[gene_col].head(10).values)
    top10_scores = df.loc[df[gene_col].isin(top10_genes), score_col].values
    rest_scores = df.loc[~df[gene_col].isin(top10_genes), score_col].values

    neural_mask = df["proof_status"] == "known_rnai_validated"
    neural_scores = df.loc[neural_mask, score_col].values
    non_neural_scores = df.loc[~neural_mask, score_col].values

    results = {}

    print(f"\n--- Top-10 vs Rest ---")
    print(f"  Top-10: n={len(top10_scores)}, mean={np.mean(top10_scores):.4f}, median={np.median(top10_scores):.4f}")
    print(f"  Rest:   n={len(rest_scores)}, mean={np.mean(rest_scores):.4f}, median={np.median(rest_scores):.4f}")

    cd_top10 = cliffs_delta(top10_scores, rest_scores)
    d_top10 = cohens_d(top10_scores, rest_scores)
    g_top10 = hedges_g(top10_scores, rest_scores)
    u_top10, p_top10 = stats.mannwhitneyu(top10_scores, rest_scores, alternative="greater")

    results["top10_vs_rest"] = {
        "caveat": "tautological: groups defined by the tested score",
        "cliffs_delta": float(cd_top10),
        "cohens_d": float(d_top10),
        "hedges_g": float(g_top10),
        "mann_whitney_u": float(u_top10),
        "p_value": float(p_top10),
        "top10_n": int(len(top10_scores)),
        "rest_n": int(len(rest_scores)),
        "top10_mean": float(np.mean(top10_scores)),
        "rest_mean": float(np.mean(rest_scores)),
    }
    print(f"  Cliff's delta: {cd_top10:.4f}")
    print(f"  Cohen's d:     {d_top10:.4f}")
    print(f"  Hedges' g:     {g_top10:.4f}")
    print(f"  Mann-Whitney U: {u_top10:.1f}, p={p_top10:.4e}")

    print(f"\n--- Neural vs Non-Neural ---")
    print(f"  Neural:     n={len(neural_scores)}, mean={np.mean(neural_scores):.4f}")
    print(f"  Non-neural: n={len(non_neural_scores)}, mean={np.mean(non_neural_scores):.4f}")
    print(f"  (circular: the rnai stream IS the proof_status label)")

    cd_neural = cliffs_delta(neural_scores, non_neural_scores)
    d_neural = cohens_d(neural_scores, non_neural_scores)
    g_neural = hedges_g(neural_scores, non_neural_scores)
    u_neural, p_neural = stats.mannwhitneyu(neural_scores, non_neural_scores, alternative="greater")

    results["neural_vs_non_neural"] = {
        "caveat": "circular: the rnai stream encodes the proof_status label",
        "cliffs_delta": float(cd_neural),
        "cohens_d": float(d_neural),
        "hedges_g": float(g_neural),
        "mann_whitney_u": float(u_neural),
        "p_value": float(p_neural),
        "neural_n": int(len(neural_scores)),
        "non_neural_n": int(len(non_neural_scores)),
        "neural_mean": float(np.mean(neural_scores)),
        "non_neural_mean": float(np.mean(non_neural_scores)),
    }
    print(f"  Cliff's delta: {cd_neural:.4f}")
    print(f"  Cohen's d:     {d_neural:.4f}")
    print(f"  Hedges' g:     {g_neural:.4f}")
    print(f"  Mann-Whitney U: {u_neural:.1f}, p={p_neural:.4e}")

    # Honest contrast: same neural label, score WITHOUT label-derived
    # streams (rnai, neural_enriched, neural_specificity).
    hs = honest_score(df)
    neural_h = hs[neural_mask].values
    non_neural_h = hs[~neural_mask].values
    cd_h = cliffs_delta(neural_h, non_neural_h)
    d_h = cohens_d(neural_h, non_neural_h)
    g_h = hedges_g(neural_h, non_neural_h)
    u_h, p_h = stats.mannwhitneyu(neural_h, non_neural_h, alternative="greater")
    results["neural_vs_non_neural_honest"] = {
        "caveat": ("honest: score recomputed excluding rnai/neural_* "
                   "label-derived streams"),
        "cliffs_delta": float(cd_h),
        "cohens_d": float(d_h),
        "hedges_g": float(g_h),
        "mann_whitney_u": float(u_h),
        "p_value": float(p_h),
        "neural_n": int(len(neural_h)),
        "non_neural_n": int(len(non_neural_h)),
        "neural_mean": float(np.mean(neural_h)),
        "non_neural_mean": float(np.mean(non_neural_h)),
    }
    print(f"\n--- Neural vs Non-Neural (honest score) ---")
    print(f"  Neural: mean={np.mean(neural_h):.4f}  Non-neural: mean={np.mean(non_neural_h):.4f}")
    print(f"  Cliff's delta: {cd_h:.4f}")
    print(f"  Cohen's d:     {d_h:.4f}")
    print(f"  Hedges' g:     {g_h:.4f}")
    print(f"  Mann-Whitney U: {u_h:.1f}, p={p_h:.4e}")

    out_path = RESULTS_DIR / "effect_sizes.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
