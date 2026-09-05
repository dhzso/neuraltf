#!/usr/bin/env python
"""Mann-Whitney U test for top-10 vs remaining candidates.

2026-09-04 audit note (selection tautology, documented): the "top-10"
group is defined as the 10 highest values of the SAME score being
tested, so the U-test on score vs selection-by-score is tautological -
under any score distribution the top-10 exceed the rest, U = n1*n2 and
delta = 1 by construction. The U/rank-biserial/Cliff's-delta FORMULAS
are reported for completeness but carry no evidential value for the
selection itself.

The script therefore reports two comparisons:
1. circular: top-10 by integrated_score vs rest (tautological; kept for
   continuity with figure 23's dual circular/honest framing).
2. honest: the top-10 shortlist vs rest on a RECOMPUTED score that
   EXCLUDES the streams directly derived from the ground-truth label
   (rnai, neural_enriched, neural_specificity) - the label-leakage-free
   contrast (same recomputation as precision_recall.py).

Usage:
    python scripts/stats/mann_whitney_top10.py
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
# Label-derived streams excluded from the honest score (they encode the
# RNAi/neural ground truth the shortlist was partly selected on):
HONEST_EXCLUDE = ["rnai", "neural_enriched", "neural_specificity"]


def rank_biserial_correlation(u_stat, n1, n2):
    """Rank-biserial correlation from Mann-Whitney U.

    U counts pairs where a group-1 value EXCEEDS a group-2 value (plus
    half the ties). The conventional effect size favoring group 1 is
    r = 2U/(n1*n2) - 1, giving +1 for perfect separation (group 1 all
    higher) and -1 for the reverse.
    """
    return (2.0 * u_stat) / (n1 * n2) - 1.0


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


def compare(top_scores, rest_scores):
    n1, n2 = len(top_scores), len(rest_scores)
    u_stat, p_one = stats.mannwhitneyu(top_scores, rest_scores, alternative="greater")
    _, p_two = stats.mannwhitneyu(top_scores, rest_scores, alternative="two-sided")
    r_rb = rank_biserial_correlation(u_stat, n1, n2)
    cliff = 0.0
    for t in top_scores:
        cliff += np.sum(t > rest_scores) - np.sum(t < rest_scores)
    cliff /= (n1 * n2)
    return {
        "top10_n": int(n1), "rest_n": int(n2),
        "top10_mean": float(np.mean(top_scores)),
        "rest_mean": float(np.mean(rest_scores)),
        "top10_median": float(np.median(top_scores)),
        "rest_median": float(np.median(rest_scores)),
        "u_statistic": float(u_stat),
        "p_value_one_sided": float(p_one),
        "p_value_two_sided": float(p_two),
        "rank_biserial_correlation": float(r_rb),
        "cliffs_delta": float(cliff),
    }


def main():
    print("=== Mann-Whitney U: Top-10 vs Rest (circular + honest) ===")

    p = RUN_DIR / "rank.csv"
    if not p.exists():
        print("Error: rank.csv not found")
        return 1
    df = pd.read_csv(p).drop_duplicates(subset="gene_id", keep="first")

    score_col = "integrated_score"
    gene_col = "gene_id"
    df = df.dropna(subset=[score_col])

    df_sorted = df.sort_values(score_col, ascending=False)
    top10_genes = set(df_sorted[gene_col].head(10).values)

    # 1) circular comparison (tautological by construction)
    top10_scores = df.loc[df[gene_col].isin(top10_genes), score_col].values
    rest_scores = df.loc[~df[gene_col].isin(top10_genes), score_col].values
    circular = compare(top10_scores, rest_scores)
    print(f"\n[circular] Top-10: n={circular['top10_n']}, "
          f"mean={circular['top10_mean']:.4f}  |  Rest: n={circular['rest_n']}, "
          f"mean={circular['rest_mean']:.4f}")
    print(f"[circular] U={circular['u_statistic']:.0f}  p={circular['p_value_one_sided']:.3e}  "
          f"r_rb={circular['rank_biserial_correlation']:.3f}  delta={circular['cliffs_delta']:.3f}"
          f"  (tautological: groups defined by the tested score)")

    # 2) honest comparison: same shortlist, label-independent score
    hs = honest_score(df)
    top10_h = hs[df[gene_col].isin(top10_genes)].values
    rest_h = hs[~df[gene_col].isin(top10_genes)].values
    honest = compare(top10_h, rest_h)
    print(f"\n[honest]   Top-10: mean={honest['top10_mean']:.4f}  |  "
          f"Rest: mean={honest['rest_mean']:.4f}")
    print(f"[honest]   U={honest['u_statistic']:.0f}  p={honest['p_value_one_sided']:.3e}  "
          f"r_rb={honest['rank_biserial_correlation']:.3f}  delta={honest['cliffs_delta']:.3f}"
          f"  (score excludes rnai/neural_* label-derived streams)")

    results = {
        "test": "Mann-Whitney U (one-sided, greater)",
        "caveat": (
            "The 'circular' comparison selects the top-10 by the tested "
            "score itself - tautological (U=n1*n2, delta=1 by "
            "construction). The 'honest' comparison recomputes the score "
            "excluding rnai/neural_enriched/neural_specificity "
            "(label-derived streams)."
        ),
        "circular": circular,
        "honest": honest,
        "top10_genes": sorted(list(top10_genes)),
    }

    print(f"\nTop-10 genes (by integrated_score):")
    for _, row in df_sorted.head(10).iterrows():
        nm = row.get("gene_name", row[gene_col])
        nm = nm if isinstance(nm, str) and str(nm) != "nan" else row[gene_col]
        print(f"  {str(nm)[:28]:>28}  score={row[score_col]:.4f}")

    out_path = RESULTS_DIR / "mann_whitney_top10.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

