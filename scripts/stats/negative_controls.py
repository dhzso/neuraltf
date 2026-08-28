#!/usr/bin/env python
"""Negative control analysis for scoring specificity.

Scores 100 random non-TF genes and 100 random non-neural TFs as negative
controls. Compares their score distributions to neural TF candidates.

Usage:
    python scripts/stats/negative_controls.py --n-controls 100 --seed 42
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
FIG_DIR = REPO / "projects" / "NeuralTF" / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


def load_all_genes():
    """Load all candidates to identify TFs and non-TFs."""
    for fname in ["supplementary_table_S2_fixed_all249.csv", "fstf_ranked_all.csv"]:
        p = RESULTS_DIR / fname
        if p.exists():
            return pd.read_csv(p)
    p = RUN_DIR / "rank.csv"
    if p.exists():
        return pd.read_csv(p)
    raise FileNotFoundError("No candidate score file found")


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

    neural_mask = df["proof_status"] == "known_rnai_validated"
    tf_mask = df.get("mmc4_tf_flag", df.get("perez_tf_class", pd.Series(dtype=str))).notna()

    neural_scores = df.loc[neural_mask, score_col].dropna().values
    all_non_neural = df.loc[~neural_mask, score_col].dropna().values

    tf_non_neural = df.loc[~neural_mask & tf_mask, score_col].dropna().values
    non_tf_scores = df.loc[~neural_mask & ~tf_mask, score_col].dropna().values

    rng = np.random.default_rng(args.seed)

    n_ctrl = min(args.n_controls, len(all_non_neural))
    random_non_tf_idx = rng.choice(len(non_tf_scores), size=min(n_ctrl, len(non_tf_scores)), replace=False) if len(non_tf_scores) > 0 else np.array([], dtype=int)
    random_non_neural_tf_idx = rng.choice(len(tf_non_neural), size=min(n_ctrl, len(tf_non_neural)), replace=False) if len(tf_non_neural) > 0 else np.array([], dtype=int)

    ctrl_non_tf = non_tf_scores[random_non_tf_idx] if len(random_non_tf_idx) > 0 else np.array([])
    ctrl_non_neural_tf = tf_non_neural[random_non_neural_tf_idx] if len(random_non_neural_tf_idx) > 0 else np.array([])

    print(f"Neural TFs (ground truth): {len(neural_scores)}")
    print(f"Non-TF controls: {len(ctrl_non_tf)}")
    print(f"Non-neural TF controls: {len(ctrl_non_neural_tf)}")

    results = {}
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

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    data_plot = [neural_scores, ctrl_non_tf] if len(ctrl_non_tf) > 0 else [neural_scores]
    labels_plot = ["Neural TFs", "Random Non-TFs"] if len(ctrl_non_tf) > 0 else ["Neural TFs"]
    bp = ax.boxplot(data_plot, labels=labels_plot, patch_artist=True,
                    boxprops=dict(facecolor="#4C72B0", alpha=0.7),
                    medianprops=dict(color="red"))
    if len(data_plot) > 1:
        bp["boxes"][1].set_facecolor("#55A868")
    ax.set_ylabel("Integrated Score")
    ax.set_title("Neural TFs vs Random Non-TFs")
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    data_plot2 = [neural_scores, ctrl_non_neural_tf] if len(ctrl_non_neural_tf) > 0 else [neural_scores]
    labels_plot2 = ["Neural TFs", "Non-Neural TFs"] if len(ctrl_non_neural_tf) > 0 else ["Neural TFs"]
    bp2 = ax.boxplot(data_plot2, labels=labels_plot2, patch_artist=True,
                     boxprops=dict(facecolor="#4C72B0", alpha=0.7),
                     medianprops=dict(color="red"))
    if len(data_plot2) > 1:
        bp2["boxes"][1].set_facecolor("#C44E52")
    ax.set_ylabel("Integrated Score")
    ax.set_title("Neural TFs vs Non-Neural TFs")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = FIG_DIR / "negative_controls.png"
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)
    print(f"Saved: {fig_path}")

    out_path = RESULTS_DIR / "negative_control_stats.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
