#!/usr/bin/env python
"""Negative control analysis for scoring specificity.

2026-09-06 audit redesign (the 2026-09-04 availability-matching pass was
necessary but not sufficient):

1. LABEL-FREE TEST SCORE. The previous headline test ran on
   integrated_score, which contains the `rnai` stream — a perfect 1/0
   copy of the group label (rnai==1.0 for ALL validated genes, 0.0 for
   everyone else). Any contrast on that score is guaranteed-separable
   by construction. The test score now EXCLUDES rnai, neural_enriched,
   neural_specificity, and perez_lineage (a hand-curated neural-family
   list with AUC 0.99 alone — the same species of leakage as rnai).
2. HONEST TF CONTROLS. Availability-matching on n_streams==7 dragged
   neural_enriched==1 genes into the "non-neural TF" control group
   (48/100 in the shipped run — the only route past the 5 always-present
   streams is the neural_* pair). The strict control pool now
   additionally requires neural_enriched != 1, and the script reports
   the within-neural contamination fraction of the permissive pool so
   neither number can be silently over-read.
3. LOUD FALLBACKS. Empty control pools previously fell back to
   bottom-quartile-by-tested-score genes (selection on the outcome) or
   all non-neural genes, silently. Both now hard-fail with a message.
4. Draw uncertainty. The matched draw is repeated over n_draws seeds;
   the reported p is the median across draws with a dispersion
   quantile, so a single unlucky/lucky draw can no longer decide the
   headline.

Groups:
  - "neural":    RNAi-screened (proof_status) — King mmc5 screening-list
                 genes (phenotype status tracked separately).
  - "non_tf":    candidates with no Perez TF-class evidence
                 (perez_lineage == 0/NaN) — lowest-confidence controls.
  - "non_neural_tf": TF-classified candidates WITHOUT the screened
                 label AND without King-neural enrichment — the
                 strictest like-for-like control.

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
           "perez_lineage", "perez_influence", "fincher_brain", "cui_temporal"]
W_DEFAULT = np.array([0.1, 0.1, 0.1, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
# Streams that encode the RNAi/neural labels (or curated neural-family
# membership). NEVER allowed into the tested score.
LEAKING_STREAMS = ["rnai", "neural_enriched", "neural_specificity",
                   "perez_lineage"]


def label_free_score(df: pd.DataFrame) -> pd.Series:
    """Renormalized weighted score over label-independent streams only
    (identical to EvidenceScorer's per-record renormalization)."""
    keep = [s for s in STREAMS if s in df.columns and s not in LEAKING_STREAMS]
    W = np.array([W_DEFAULT[STREAMS.index(s)] for s in keep])
    M = df[keep].to_numpy(dtype=float)
    mask = ~np.isnan(M)
    num = np.where(mask, M, 0.0) @ W
    den = (mask * W).sum(axis=1)
    return pd.Series(np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0),
                     index=df.index)


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

    Genes are ranked by |n_streams_avail - target| (ties broken by a random
    key so the draw stays random within the matched band) and the top
    n_draw are taken. This removes the availability confound: validated
    TFs have ~9 non-null streams; controls with ~5 would trivially
    separate on any score that renormalizes over present streams.
    """
    if pool.empty:
        return pool
    pool = pool.copy()
    match_col = "n_streams_avail" if "n_streams_avail" in pool.columns else "n_streams"
    pool["_dist"] = (pool[match_col] - target_n_streams).abs()
    pool["_rand"] = rng.random(len(pool))
    pool = pool.sort_values(["_dist", "_rand"])
    take = pool.head(min(n_draw, len(pool)))
    return take.drop(columns=["_dist", "_rand"])


def main():
    parser = argparse.ArgumentParser(description="Negative control analysis")
    parser.add_argument("--n-controls", type=int, default=100, help="Number of random controls per group")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (base; draw-uncertainty uses seed..seed+n_draws-1)")
    parser.add_argument("--n-draws", type=int, default=20, help="Matched-draw repetitions for draw-uncertainty reporting")
    args = parser.parse_args()

    print("=== Negative Control Analysis (label-free score, matched, draw-uncertainty) ===")

    df = load_all_genes()

    # TEST SCORE: label-free (excludes rnai/neural_*/perez_lineage).
    # The tested score must NEVER contain a copy of the group label.
    df["_score"] = label_free_score(df)
    score_col = "_score"

    # stream availability per gene (matching covariate)
    stream_cols = [s for s in STREAMS if s in df.columns]
    df["n_streams_avail"] = df[stream_cols].notna().sum(axis=1)

    gene_col = "gene_id" if "gene_id" in df.columns else "gene_id_v6"

    # Neural candidates: RNAi-screened ONLY (neural_enriched is a scoring
    # stream — including it would make the control circular)
    neural_mask = df["proof_status"] == "tested"
    # TF vs non-TF indicator
    if "perez_lineage" in df.columns:
        tf_mask = df["perez_lineage"].fillna(0) > 0
    else:
        tf_mask = df["proof_status"].notna()

    neural_df = df.loc[neural_mask].dropna(subset=[score_col])
    target_avail = int(np.median(neural_df["n_streams_avail"])) if len(neural_df) else 9

    pool_non_tf = df.loc[~neural_mask & ~tf_mask].dropna(subset=[score_col])
    # STRICT TF control pool: also excludes King-neural-enriched genes.
    # Availability-matching alone dragged neural_enriched==1 genes into
    # this pool (48/100 in the 2026-09-04 run) because the neural_* pair
    # is the only route past the 5 always-present streams — the
    # "non-neural TF" contrast was half within-neural.
    king_neural = (df["neural_enriched"].fillna(0) > 0) if "neural_enriched" in df.columns \
        else pd.Series(False, index=df.index)
    pool_tf_non_neural = df.loc[~neural_mask & tf_mask & ~king_neural].dropna(subset=[score_col])
    # Permissive pool (neural-enriched TFs allowed) retained only to
    # REPORT the contamination fraction, never as the tested control.
    pool_tf_permissive = df.loc[~neural_mask & tf_mask].dropna(subset=[score_col])
    if len(pool_tf_permissive):
        contam = float((pool_tf_permissive["neural_enriched"].fillna(0) > 0).mean()) \
            if "neural_enriched" in pool_tf_permissive.columns else float("nan")
    else:
        contam = float("nan")

    # LOUD fallbacks (the previous silent fallbacks selected controls on
    # the tested score — guaranteed separation — or widened to all
    # non-neural genes without any warning).
    if pool_non_tf.empty:
        print("ERROR: no non-TF control pool exists (every non-validated "
              "candidate carries a Perez TF class). Refusing the previous "
              "silent bottom-quartile-by-score fallback (it selects on "
              "the tested score). Provide a genuine control pool.")
        return 1
    if pool_tf_non_neural.empty:
        print("ERROR: no TF-classified, non-validated, non-neural-enriched "
              "control pool exists. Refusing the previous silent "
              "all-non-neural fallback. Provide a genuine control pool.")
        return 1

    neural_scores = neural_df[score_col].values
    if len(neural_scores) == 0:
        print("ERROR: no RNAi-screened genes with scores; nothing to test")
        return 1

    # ---- Draw-uncertainty: repeat the matched draw across seeds -------
    p1_draws, d1_draws = [], []
    p2_draws, d2_draws = [], []
    ctrl_non_tf_keep = ctrl_tf_keep = None
    for draw in range(args.n_draws):
        rng = np.random.default_rng(args.seed + draw)
        ctrl_non_tf_df = matched_sample(pool_non_tf, target_avail,
                                        args.n_controls, rng)
        ctrl_tf_df = matched_sample(pool_tf_non_neural, target_avail,
                                    args.n_controls, rng)
        if draw == 0:
            ctrl_non_tf_keep = ctrl_non_tf_df
            ctrl_tf_keep = ctrl_tf_df
        c1 = ctrl_non_tf_df[score_col].values
        c2 = ctrl_tf_df[score_col].values
        if len(c1):
            _, p1 = stats.mannwhitneyu(neural_scores, c1, alternative="greater")
            p1_draws.append(float(p1))
            d1_draws.append(pooled_cohens_d(neural_scores, c1))
        if len(c2):
            _, p2 = stats.mannwhitneyu(neural_scores, c2, alternative="greater")
            p2_draws.append(float(p2))
            d2_draws.append(pooled_cohens_d(neural_scores, c2))

    ctrl_non_tf_df = ctrl_non_tf_keep
    ctrl_tf_df = ctrl_tf_keep
    ctrl_non_tf = ctrl_non_tf_df[score_col].values if ctrl_non_tf_df is not None and len(ctrl_non_tf_df) else np.array([])
    ctrl_non_neural_tf = ctrl_tf_df[score_col].values if ctrl_tf_df is not None and len(ctrl_tf_df) else np.array([])

    print(f"Test score: label-free (excludes {'/'.join(LEAKING_STREAMS)})")
    print(f"Neural TFs (King mmc5 RNAi-screened label): {len(neural_scores)} "
          f"(median streams available: {target_avail})")
    print(f"Non-TF controls (matched, no Perez TF class): {len(ctrl_non_tf)} "
          f"(median streams: "
          f"{int(np.median(ctrl_non_tf_df['n_streams_avail'])) if len(ctrl_non_tf_df) else '-'})")
    print(f"Non-neural TF controls (matched, TF-classified, neural_enriched!=1): "
          f"{len(ctrl_non_neural_tf)} "
          f"(median streams: "
          f"{int(np.median(ctrl_tf_df['n_streams_avail'])) if len(ctrl_tf_df) else '-'})")
    print(f"Permissive TF pool neural-enriched contamination (reported, not tested): "
          f"{contam:.1%}")

    results = {}

    def _draw_summary(p_draws, d_draws):
        if not p_draws:
            return None
        return {
            "p_median": float(np.median(p_draws)),
            "p_q25": float(np.percentile(p_draws, 25)),
            "p_q75": float(np.percentile(p_draws, 75)),
            "p_max_draw_p": float(max(p_draws)),
            "d_median": float(np.median(d_draws)) if d_draws else None,
            "n_draws": len(p_draws),
        }

    if len(ctrl_non_tf) > 0:
        u1, p1 = stats.mannwhitneyu(neural_scores, ctrl_non_tf, alternative="greater")
        d1 = pooled_cohens_d(neural_scores, ctrl_non_tf)
        results["neural_vs_random_non_tf"] = {
            "mann_whitney_u": float(u1),
            "p_value": float(p1),
            "cohens_d": float(d1),
            "neural_mean": float(np.mean(neural_scores)),
            "control_mean": float(np.mean(ctrl_non_tf)),
            "matched_on": "n_streams_avail",
            "score": "label-free (rnai/neural_*/perez_lineage excluded)",
            "draw_uncertainty": _draw_summary(p1_draws, d1_draws),
        }
        print(f"  Neural vs Random Non-TF: U={u1:.1f}, p={p1:.4e}, d={d1:.3f} "
              f"(median p across {len(p1_draws)} draws: {np.median(p1_draws):.4e})")

    if len(ctrl_non_neural_tf) > 0:
        u2, p2 = stats.mannwhitneyu(neural_scores, ctrl_non_neural_tf, alternative="greater")
        d2 = pooled_cohens_d(neural_scores, ctrl_non_neural_tf)
        results["neural_vs_non_neural_tf"] = {
            "mann_whitney_u": float(u2),
            "p_value": float(p2),
            "cohens_d": float(d2),
            "neural_mean": float(np.mean(neural_scores)),
            "control_mean": float(np.mean(ctrl_non_neural_tf)),
            "matched_on": "n_streams_avail; neural_enriched!=1",
            "score": "label-free (rnai/neural_*/perez_lineage excluded)",
            "draw_uncertainty": _draw_summary(p2_draws, d2_draws),
        }
        print(f"  Neural vs Non-Neural TF: U={u2:.1f}, p={p2:.4e}, d={d2:.3f} "
              f"(median p across {len(p2_draws)} draws: {np.median(p2_draws):.4e})")

    # JSON export with keys needed by figure 24. Empty groups export an
    # EMPTY list (never the [0.0] sentinel — a phantom 0-score gene).
    results["neural_tfs"] = [float(x) for x in neural_scores[:100]] if len(neural_scores) else []
    results["non_tfs"] = [float(x) for x in ctrl_non_tf[:100]] if len(ctrl_non_tf) else []
    # key kept as "random" for figure-24 compatibility, but these are
    # availability-matched, neural-free TF controls (not permutations)
    results["random"] = [float(x) for x in ctrl_non_neural_tf[:100]] if len(ctrl_non_neural_tf) else []
    results["group_labels"] = {
        "neural_tfs": "RNAi-screened TFs (King mmc5 screening list; phenotype status tracked separately)",
        "non_tfs": "availability-matched candidates without Perez TF class",
        "random": "availability-matched TF-classified, non-validated, non-neural-enriched candidates",
    }
    results["p_neural_vs_non"] = float(results.get("neural_vs_random_non_tf", {}).get("p_value", 1.0))
    results["score_note"] = (
        "All contrasts run on the LABEL-FREE score (rnai, neural_enriched, "
        "neural_specificity, perez_lineage excluded). The rnai stream is a "
        "perfect 1/0 copy of the group label — the previous headline "
        "p-values (e.g. 1.14e-25) were partly arithmetically guaranteed."
    )
    results["permissive_pool_neural_contamination"] = float(contam) if contam == contam else None

    # NOTE: no standalone PNG — the numbered publication figure
    # (figures/24_negative_controls.py) renders the showcase version.

    # ---- 2026-09-11: phenotype-confirmed label arm ----------------------
    # 'tested' = King mmc5 SCREENED list ("All Transcription Factors
    # Inhibited"; the distributed copy lost the red/green phenotype font
    # encoding). The FISH-confirmed subset (paper Fig 3J/4E, S4, S7, S8)
    # is the strictest defensible positive group; rerun the primary
    # contrast against it with the same matched controls.
    from bioforge.evidence.groundtruth import PHENOTYPE_CONFIRMED_V6
    conf_mask = (neural_mask
                 & df[gene_col].astype(str).isin(PHENOTYPE_CONFIRMED_V6)
                 & df[score_col].notna())
    conf_scores = df.loc[conf_mask, score_col].values
    results["phenotype_confirmed_arm"] = {
        "n_positives": int(len(conf_scores)),
        "label_note": ("FISH-confirmed loss-of-cell-type phenotypes "
                       "(King 2024 Fig 3J/4E, S4, S7, S8)"),
    }
    if len(conf_scores) >= 2 and len(ctrl_non_tf) > 0:
        u3, p3 = stats.mannwhitneyu(conf_scores, ctrl_non_tf, alternative="greater")
        d3 = pooled_cohens_d(conf_scores, ctrl_non_tf)
        results["phenotype_confirmed_arm"]["vs_non_tf"] = {
            "mann_whitney_u": float(u3),
            "p_value": float(p3),
            "cohens_d": float(d3),
            "positive_mean": float(np.mean(conf_scores)),
            "control_mean": float(np.mean(ctrl_non_tf)),
        }
        print(f"  Phenotype-confirmed vs Random Non-TF: U={u3:.1f}, "
              f"p={p3:.4e}, d={d3:.3f}")
    if len(conf_scores) >= 2 and len(ctrl_non_neural_tf) > 0:
        u4, p4 = stats.mannwhitneyu(conf_scores, ctrl_non_neural_tf,
                                    alternative="greater")
        d4 = pooled_cohens_d(conf_scores, ctrl_non_neural_tf)
        results["phenotype_confirmed_arm"]["vs_non_neural_tf"] = {
            "mann_whitney_u": float(u4),
            "p_value": float(p4),
            "cohens_d": float(d4),
            "positive_mean": float(np.mean(conf_scores)),
            "control_mean": float(np.mean(ctrl_non_neural_tf)),
        }
        print(f"  Phenotype-confirmed vs Non-Neural TF: U={u4:.1f}, "
              f"p={p4:.4e}, d={d4:.3f}")

    results["label_note"] = (
        "2026-09-11: the 'neural' group (proof_status=='tested') is the "
        "King mmc5 RNAi SCREENING list — screened, not phenotype-validated. "
        "The phenotype_confirmed_arm keys hold the stricter-label contrasts."
    )

    out_path = RESULTS_DIR / "negative_control_stats.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
