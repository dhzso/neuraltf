#!/usr/bin/env python
"""Score-shuffling permutation test (alternative null model).

2026-09-06 audit redesign (two broken aspects of the previous null):

1. JOINT whole-row permutation. The previous null shuffled each stream
   column INDEPENDENTLY. But within a gene the streams are determin-
   istically dependent by pipeline construction (reproducibility is a
   function of atlas membership which co-occurs with expression/
   specificity; neural_enriched/neural_specificity are always jointly
   present; perez_influence availability is gated on Perez membership).
   Independent per-column shuffles destroy this dependence, so real
   genes — even under "no biology" — get correlated stream values the
   null can never produce. The permutation therefore partly tested
   "streams co-occur within genes" (true by construction), inflating
   significance (~1,200/11,675 at p<0.05 against ~584 expected).
   The null now permutes WHOLE ROWS (stream vectors) across candidates,
   preserving the within-gene dependence structure exactly: H0 = "the
   association between genes and their evidence vectors is exchangeable".
2. Family multiplicity. 11,675 per-gene p-values are now BH-corrected
   across the full family, and the resolution floor 1/(n+1) is stated
   against the family threshold (per-gene alpha was never the operative
   level; nothing downstream may quote uncorrected p's).

Usage:
    python scripts/stats/score_shuffling_permutation.py --n-perm 1000 --seed 42
"""

import argparse
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

STREAMS = [
    "expression", "specificity", "reproducibility", "rnai",
    "correlation", "neural_enriched", "neural_specificity",
    "perez_lineage", "perez_influence",
]
# Must match bioforge.evidence.scoring.DEFAULT_WEIGHTS exactly
# (the old perez_lineage=0.2 made "real" scores differ from rank.csv).
W_DEFAULT = np.array([0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])


def compute_all_integrated_scores(scores: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Vectorized calculation of integrated scores across all candidates."""
    scores_filled = np.nan_to_num(scores, nan=0.0)
    valid_mask = ~np.isnan(scores)
    w_denom = np.dot(valid_mask, weights)
    w_denom = np.where(w_denom > 0, w_denom, 1.0)
    return np.dot(scores_filled, weights) / w_denom


def shuffle_rows(scores_matrix, rng):
    """Permute WHOLE evidence-vector rows across candidates (joint
    permutation preserving within-gene stream dependence).

    H0: the assignment of complete evidence vectors to genes is
    exchangeable. Permuting rows (not columns) keeps every within-gene
    co-occurrence pattern intact — neural_enriched stays with its
    neural_specificity, reproducibility stays with the expression/
    specificity profile it was derived from — so the null isolates the
    gene-evidence ASSOCIATION instead of also testing deterministic
    within-gene construction (which the previous per-column shuffles
    violated).
    """
    return scores_matrix[rng.permutation(scores_matrix.shape[0])]


def main():
    parser = argparse.ArgumentParser(description="Score-shuffling permutation test")
    parser.add_argument("--n-perm", type=int, default=1000, help="Number of permutations (default: 1000)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--candidates", action="store_true",
                        help="restrict p-value aggregation and family BH to neural candidates (rank_neural.csv) — the targeted high-resolution mode")
    args = parser.parse_args()

    print(f"=== Score-Shuffling Permutation Test (n={args.n_perm}{', --candidates mode' if args.candidates else ''}) ===")

    rng = np.random.default_rng(args.seed)

    p = RUN_DIR / "rank.csv"
    if not p.exists():
        print("Error: no candidate score file found (run the pipeline first)")
        return 1

    df = pd.read_csv(p).drop_duplicates(subset="gene_id", keep="first")

    gene_col = "gene_id" if "gene_id" in df.columns else "gene_id_v6"
    score_col = "integrated_score" if "integrated_score" in df.columns else "composite_score"

    stream_cols = [s for s in STREAMS if s in df.columns]
    if len(stream_cols) < 3:
        stream_cols = [c for c in df.columns if c in STREAMS]

    scores_matrix = df[stream_cols].values.astype(float)
    weights = W_DEFAULT[:len(stream_cols)]
    weights = weights / weights.sum()

    real_scores = compute_all_integrated_scores(scores_matrix, weights)
    df["real_integrated_score"] = real_scores
    gene_indices = {gene: idx for idx, gene in enumerate(df[gene_col].values)}

    if args.candidates:
        neural_csv = RUN_DIR / "rank_neural.csv"
        if not neural_csv.exists():
            print(f"--candidates requires {neural_csv}; run pipeline first")
            return 1
        neural_ids = set(pd.read_csv(neural_csv)["gene_id"].astype(str))
        target_genes = [g for g in df[gene_col].values if g in neural_ids]
        print(f"--candidates mode: restricting p-value family to {len(target_genes)} neural candidates")
    else:
        target_genes = list(df[gene_col].values)

    df_target = df[df[gene_col].isin(set(target_genes))].sort_values("real_integrated_score", ascending=False)

    print(f"Universe: {len(df)} genes, Target family: {len(target_genes)} genes, Streams: {len(stream_cols)}")
    print(f"Top-10 real scores in target: {df_target['real_integrated_score'].head(10).values}")

    null_distributions = {gene: [] for gene in target_genes}
    target_idx = [gene_indices[g] for g in target_genes]
    family_real_mean = float(np.mean(df_target["real_integrated_score"]))
    family_null_means = []

    for perm in range(args.n_perm):
        shuffled = shuffle_rows(scores_matrix, rng)
        perm_scores = compute_all_integrated_scores(shuffled, weights)
        for g, idx in zip(target_genes, target_idx):
            null_distributions[g].append(perm_scores[idx])
        family_null_means.append(float(np.mean([perm_scores[idx] for idx in target_idx])))

        if (perm + 1) % 200 == 0 or (perm + 1) == args.n_perm:
            print(f"  Completed {perm+1}/{args.n_perm} permutations")

    print("\n=== Computing Empirical P-values (joint row permutation) ===")
    pvals = []
    for _, row in df_target.iterrows():
        gene = row[gene_col]
        real_s = row["real_integrated_score"]
        null_dist = np.array(null_distributions[gene])
        p = (np.sum(null_dist >= real_s) + 1) / (args.n_perm + 1)
        pvals.append(p)

    # BH-FDR across the target family (134 genes in --candidates mode, 11,675 in full mode)
    from statsmodels.stats.multitest import multipletests
    _, qvals, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")

    df_out = df_target[[gene_col, "gene_name", "real_integrated_score"]].copy() if "gene_name" in df_target.columns else df_target[[gene_col, "real_integrated_score"]].copy()
    df_out = df_out.reset_index(drop=True)
    df_out["empirical_p_shuffled"] = pvals
    df_out["q_bh_family"] = qvals
    df_out["significant_fdr_05"] = qvals < 0.05

    # Non-degenerate family-mean statistic: test whether the candidate family's
    # average score exceeds a random draw of equal size from the universe
    family_p = (np.sum(np.array(family_null_means) >= family_real_mean) + 1) / (args.n_perm + 1)

    out_name = "score_shuffling_pvalues_neural.csv" if args.candidates else "score_shuffling_pvalues.csv"
    out_path = RESULTS_DIR / out_name
    df_out.to_csv(out_path, index=False)
    # Also update score_shuffling_pvalues.csv if running default
    if not args.candidates:
        df_out.to_csv(RESULTS_DIR / "score_shuffling_pvalues.csv", index=False)
    print(f"Saved: {out_path}")

    n_sig_05 = sum(p < 0.05 for p in pvals)
    n_fdr_05 = int((qvals < 0.05).sum())
    p_floor = 1.0 / (args.n_perm + 1)
    print(f"\nNominal p<0.05: {n_sig_05}/{len(pvals)} "
          f"(expected under H0: ~{0.05 * len(pvals):.0f})")
    print(f"Significant after BH-FDR q<0.05 across the {len(pvals)}-gene family: "
          f"{n_fdr_05}/{len(pvals)}")
    print(f"Resolution floor at n={args.n_perm}: min p = {p_floor:.5f} "
          f"(family Bonferroni alpha = {0.05 / len(pvals):.2e} is "
          f"{'RESOLVABLE' if p_floor < 0.05 / len(pvals) else 'NOT resolvable at this n'})")
    print(f"Family-mean statistic (candidate family mean vs null draws): real={family_real_mean:.4f}, p={family_p:.4e}")

    print("\nTop-10 by score with shuffled-null p-values (with family q):")
    for i, (_, row) in enumerate(df_out.head(10).iterrows()):
        name = row.get("gene_name", row[gene_col])
        name = str(name) if (name is not None and str(name) != "nan" and str(name) != "None") else str(row[gene_col])
        print(f"  {i+1:>2d}. {name[:24]:>24}  score={row['real_integrated_score']:.4f}  "
              f"p={row['empirical_p_shuffled']:.4f}  q={row['q_bh_family']:.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
