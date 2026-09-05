#!/usr/bin/env python
"""Negative control analysis for scoring specificity.

Groups (WS3 relabel — old labels were misleading):
  - "neural":    RNAi-validated (proof_status) — known neural TFs.
  - "non_tf":    candidates with no Perez TF-class evidence
                 (perez_lineage == 0/NaN) — lowest-confidence controls.
  - "non_neural_tf": TF-classified candidates WITHOUT the RNAi-validated
                 label — the strictest like-for-like control.

Note on circularity: `neural_enriched` is a scoring stream, so the
previous "neural" definition (validated OR neural_enriched>0) mixed the
outcome into the grouping; the RNAi-validated label alone is the honest
ground truth here.

2026-09-04 audit fix (matching): controls were previously drawn WITHOUT
matching, so group separation partly measured evidence-AVAILABILITY
differences (validated genes have more non-null streams by construction)
rather than biology. Controls are now matched to the neural group on
stream availability (nearest n_streams) before the random draw, and the
Cohen's d uses the standard pooled-SD (ddof=1) formula consistent with
effect_sizes.py (the old average-variance ddof=0 form disagreed with the
rest of the suite).

Usage:
    python scripts/stats/negative_controls.py --n-controls 100 --seed 42
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

STREAMS = ["expression", "specificity", "reproducibility", "rnai",
           "correlation", "neural_enriched", "neural_specificity",
           "perez_lineage", "perez_influence"]


def load_all_genes():
    """Load all candidates to identify TFs and non-TFs."""
    p = RUN_DIR / "rank.csv"
    if p.exists():
        return pd.read_csv(p).drop_duplicates(subset="gene_id", keep="first")
    raise FileNotFoundError("No candidate score file found in runs/pipeline_run/rank.csv")


def pooled_cohens_d(x, y):
    """Standard Cohen's d with pooled SD (ddof=1) — consistent with
    effect_sizes.py (the previous average-variance ddof=0 form was a
    different estimator under the same name)."""
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return 0.0
    pooled = np.sqrt(((nx - 1) * np.var(x, ddof=1) + (ny - 1) * np.var(y, ddof=1))
                     / (nx + ny - 2))
    if pooled == 0:
        return 0.0
    return float((np.mean(x) - np.mean(y)) / pooled)


def matched_sample(pool: pd.DataFrame, target_n_streams: int,
                   n_draw: int, rng) -> pd.DataFrame:
    """Draw controls matched on stream-availability count.

    Genes are ranked by |n_streams - target| (ties broken by a random
    key so the draw stays random within the matched band) and the top
    n_draw are taken. This removes the availability confound: validated
    TFs have ~9 non-null streams; controls with ~5 would trivially
    separate on any score that renormalizes over present streams.
    """
    if pool.empty:
        return pool
    pool = pool.copy()
    pool["_dist"] = (pool["n_streams"] - target_n_streams).abs()
    pool["_rand"] = rng.random(len(pool))
    pool = pool.sort_values(["_dist", "_rand"])
    take = pool.head(min(n_draw, len(pool)))
    return take.drop(columns=["_dist", "_rand"])


def main():
    parser = argparse.ArgumentParser(description="Negative control analysis")
    parser.add_argument("--n-controls", type=int, default=100, help="Number of random controls per group")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    print("=== Negative Control Analysis (availability-matched) ===")

    df = load_all_genes()
    score_col = None
    for c in ["integrated_score", "composite_score", "dirichlet_median_score", "fixed_weight_score"]:
        if c in df.columns:
            score_col = c
            break
    if score_col is None:
        print("Error: no score column found")
        return 1

    # stream availability per gene (matching covariate)
    stream_cols = [s for s in STREAMS if s in df.columns]
    df["n_streams_avail"] = df[stream_cols].notna().sum(axis=1)

    gene_col = "gene_id" if "gene_id" in df.columns else "gene_id_v6"

    # Neural candidates: RNAi-validated ONLY (neural_enriched is a scoring
    # stream — including it would make the control circular)
    neural_mask = df["proof_status"] == "known_rnai_validated"
    # TF vs non-TF indicator
    if "perez_lineage" in df.columns:
        tf_mask = df["perez_lineage"].fillna(0) > 0
    else:
        tf_mask = df["proof_status"].notna()

    neural_df = df.loc[neural_mask].dropna(subset=[score_col])
    target_avail = int(np.median(neural_df["n_streams_avail"])) if len(neural_df) else 9

    pool_non_tf = df.loc[~neural_mask & ~tf_mask].dropna(subset=[score_col])
    pool_tf_non_neural = df.loc[~neural_mask & tf_mask].dropna(subset=[score_col])

    # Fallback if catalog is TF-only: bottom-scoring quartile as
    # empirical non-neural baseline
    if pool_non_tf.empty:
        all_non_neural = df.loc[~neural_mask].dropna(subset=[score_col])
        cutoff = np.percentile(all_non_neural[score_col], 25)
        pool_non_tf = all_non_neural[all_non_neural[score_col] <= cutoff]
    if pool_tf_non_neural.empty:
        pool_tf_non_neural = df.loc[~neural_mask].dropna(subset=[score_col])

    rng = np.random.default_rng(args.seed)

    ctrl_non_tf_df = matched_sample(pool_non_tf, target_avail, args.n_controls, rng)
    ctrl_tf_df = matched_sample(pool_tf_non_neural, target_avail, args.n_controls, rng)

    neural_scores = neural_df[score_col].values
    ctrl_non_tf = ctrl_non_tf_df[score_col].values if len(ctrl_non_tf_df) else np.array([])
    ctrl_non_neural_tf = ctrl_tf_df[score_col].values if len(ctrl_tf_df) else np.array([])

    print(f"Neural TFs (RNAi-validated ground truth): {len(neural_scores)} "
          f"(median streams available: {target_avail})")
    print(f"Non-TF controls (matched, no Perez TF class): {len(ctrl_non_tf)} "
          f"(median streams: "
          f"{int(np.median(ctrl_non_tf_df['n_streams_avail'])) if len(ctrl_non_tf_df) else '-'})")
    print(f"Non-neural TF controls (matched, TF-classified): {len(ctrl_non_neural_tf)} "
          f"(median streams: "
          f"{int(np.median(ctrl_tf_df['n_streams_avail'])) if len(ctrl_tf_df) else '-'})")

    results = {}
    p_val_report = 1.0

    if len(ctrl_non_tf) > 0 and len(neural_scores) > 0:
        u1, p1 = stats.mannwhitneyu(neural_scores, ctrl_non_tf, alternative="greater")
        d1 = pooled_cohens_d(neural_scores, ctrl_non_tf)
        results["neural_vs_random_non_tf"] = {
            "mann_whitney_u": float(u1),
            "p_value": float(p1),
            "cohens_d": float(d1),
            "neural_mean": float(np.mean(neural_scores)),
            "control_mean": float(np.mean(ctrl_non_tf)),
            "matched_on": "n_streams_avail",
        }
        p_val_report = float(p1)
        print(f"  Neural vs Random Non-TF: U={u1:.1f}, p={p1:.4e}, d={d1:.3f}")

    if len(ctrl_non_neural_tf) > 0 and len(neural_scores) > 0:
        u2, p2 = stats.mannwhitneyu(neural_scores, ctrl_non_neural_tf, alternative="greater")
        d2 = pooled_cohens_d(neural_scores, ctrl_non_neural_tf)
        results["neural_vs_non_neural_tf"] = {
            "mann_whitney_u": float(u2),
            "p_value": float(p2),
            "cohens_d": float(d2),
            "neural_mean": float(np.mean(neural_scores)),
            "control_mean": float(np.mean(ctrl_non_neural_tf)),
            "matched_on": "n_streams_avail",
        }
        print(f"  Neural vs Non-Neural TF: U={u2:.1f}, p={p2:.4e}, d={d2:.3f}")

    # JSON export with keys needed by figure 24 (labels now honest)
    results["neural_tfs"] = [float(x) for x in neural_scores[:100]] if len(neural_scores) else [0.0]
    results["non_tfs"] = [float(x) for x in ctrl_non_tf[:100]] if len(ctrl_non_tf) else [0.0]
    # key kept as "random" for figure-24 compatibility, but these are
    # availability-matched non-neural TF candidates (not permutations)
    results["random"] = [float(x) for x in ctrl_non_neural_tf[:100]] if len(ctrl_non_neural_tf) else [0.0]
    results["group_labels"] = {
        "neural_tfs": "RNAi-validated neural TFs (ground truth)",
        "non_tfs": "availability-matched candidates without Perez TF class",
        "random": "availability-matched TF-classified, non-RNAi-validated candidates",
    }
    results["p_neural_vs_non"] = float(p_val_report)

    # NOTE: no standalone PNG — the numbered publication figure
    # (figures/24_negative_controls.py) renders the showcase version.

    out_path = RESULTS_DIR / "negative_control_stats.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
