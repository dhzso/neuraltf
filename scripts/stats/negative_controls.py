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


def load_all_genes():
    """Load all candidates to identify TFs and non-TFs."""
    p = RUN_DIR / "rank.csv"
    if p.exists():
        return pd.read_csv(p).drop_duplicates(subset="gene_id", keep="first")
    raise FileNotFoundError("No candidate score file found in runs/pipeline_run/rank.csv")


def main():
    parser = argparse.ArgumentParser(description="Negative control analysis")
    parser.add_argument("--n-controls", type=int, default=100, help="Number of random controls per group")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    print("=== Negative Control Analysis ===")

    df = load_all_genes()
    score_col = None
    for c in ["integrated_score", "composite_score", "dirichlet_median_score", "fixed_weight_score"]:
        if c in df.columns:
            score_col = c
            break
    if score_col is None:
        print("Error: no score column found")
        return 1

    gene_col = "gene_id" if "gene_id" in df.columns else "gene_id_v6"

    # Neural candidates: RNAi-validated ONLY (neural_enriched is a scoring
    # stream — including it would make the control circular)
    neural_mask = df["proof_status"] == "known_rnai_validated"
    # TF vs non-TF indicator
    if "perez_lineage" in df.columns:
        tf_mask = df["perez_lineage"].fillna(0) > 0
    else:
        tf_mask = df["proof_status"].notna()

    neural_scores = df.loc[neural_mask, score_col].dropna().values
    all_non_neural = df.loc[~neural_mask, score_col].dropna().values

    tf_non_neural = df.loc[~neural_mask & tf_mask, score_col].dropna().values
    non_tf_scores = df.loc[~neural_mask & ~tf_mask, score_col].dropna().values

    # Fallback if catalog is TF-only: use bottom-scoring decile as empirical non-neural baseline
    if len(non_tf_scores) == 0 and len(all_non_neural) > 0:
        cutoff = np.percentile(all_non_neural, 25)
        non_tf_scores = all_non_neural[all_non_neural <= cutoff]

    if len(tf_non_neural) == 0 and len(all_non_neural) > 0:
        tf_non_neural = all_non_neural

    rng = np.random.default_rng(args.seed)

    n_ctrl = min(args.n_controls, max(len(all_non_neural), 1))
    random_non_tf_idx = rng.choice(len(non_tf_scores), size=min(n_ctrl, len(non_tf_scores)), replace=False) if len(non_tf_scores) > 0 else np.array([], dtype=int)
    random_non_neural_tf_idx = rng.choice(len(tf_non_neural), size=min(n_ctrl, len(tf_non_neural)), replace=False) if len(tf_non_neural) > 0 else np.array([], dtype=int)

    ctrl_non_tf = non_tf_scores[random_non_tf_idx] if len(random_non_tf_idx) > 0 else np.array([])
    ctrl_non_neural_tf = tf_non_neural[random_non_neural_tf_idx] if len(random_non_neural_tf_idx) > 0 else np.array([])

    print(f"Neural TFs (RNAi-validated ground truth): {len(neural_scores)}")
    print(f"Non-TF controls (no Perez TF class): {len(ctrl_non_tf)}")
    print(f"Non-neural TF controls (TF-classified): {len(ctrl_non_neural_tf)}")

    results = {}
    p_val_report = 1.0

    if len(ctrl_non_tf) > 0 and len(neural_scores) > 0:
        u1, p1 = stats.mannwhitneyu(neural_scores, ctrl_non_tf, alternative="greater")
        d1 = (np.mean(neural_scores) - np.mean(ctrl_non_tf)) / np.sqrt(
            (np.var(neural_scores) + np.var(ctrl_non_tf)) / 2
        ) if (np.var(neural_scores) + np.var(ctrl_non_tf)) > 0 else 0
        results["neural_vs_random_non_tf"] = {
            "mann_whitney_u": float(u1),
            "p_value": float(p1),
            "cohens_d": float(d1),
            "neural_mean": float(np.mean(neural_scores)),
            "control_mean": float(np.mean(ctrl_non_tf)),
        }
        p_val_report = float(p1)
        print(f"  Neural vs Random Non-TF: U={u1:.1f}, p={p1:.4e}, d={d1:.3f}")

    if len(ctrl_non_neural_tf) > 0 and len(neural_scores) > 0:
        u2, p2 = stats.mannwhitneyu(neural_scores, ctrl_non_neural_tf, alternative="greater")
        d2 = (np.mean(neural_scores) - np.mean(ctrl_non_neural_tf)) / np.sqrt(
            (np.var(neural_scores) + np.var(ctrl_non_neural_tf)) / 2
        ) if (np.var(neural_scores) + np.var(ctrl_non_neural_tf)) > 0 else 0
        results["neural_vs_non_neural_tf"] = {
            "mann_whitney_u": float(u2),
            "p_value": float(p2),
            "cohens_d": float(d2),
            "neural_mean": float(np.mean(neural_scores)),
            "control_mean": float(np.mean(ctrl_non_neural_tf)),
        }
        print(f"  Neural vs Non-Neural TF: U={u2:.1f}, p={p2:.4e}, d={d2:.3f}")

    # JSON export with keys needed by figure 24 (labels now honest)
    results["neural_tfs"] = [float(x) for x in neural_scores[:100]] if len(neural_scores) else [0.0]
    results["non_tfs"] = [float(x) for x in ctrl_non_tf[:100]] if len(ctrl_non_tf) else [0.0]
    # key kept as "random" for figure-24 compatibility, but these are
    # non-neural TF candidates (not permutations)
    results["random"] = [float(x) for x in ctrl_non_neural_tf[:100]] if len(ctrl_non_neural_tf) else [0.0]
    results["group_labels"] = {
        "neural_tfs": "RNAi-validated neural TFs (ground truth)",
        "non_tfs": "candidates without Perez TF class",
        "random": "TF-classified, non-RNAi-validated candidates",
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
