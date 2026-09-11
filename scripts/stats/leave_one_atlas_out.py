#!/usr/bin/env python
"""Leave-one-STREAM-out sensitivity analysis.

2026-09-04 audit relabel: despite its historical filename, this analysis
leaves out one EVIDENCE STREAM (of the 9), not one atlas - the atlas
contributions are already collapsed into single columns (expression is
max over Fincher/Plass/Cui/King), so individual-atlas attribution cannot
be isolated from rank.csv. The filename is kept for driver/figure
compatibility; the docstring, output columns, and figure labels now say
"stream".

Re-ranks candidates with each evidence stream removed (weights
renormalized over the remaining streams exactly as EvidenceScorer does),
tracking top-10 stability across leave-one-out iterations.

Usage:
    python scripts/stats/leave_one_atlas_out.py
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
    "perez_lineage", "perez_influence", "fincher_brain", "cui_temporal",
]
W_DEFAULT = np.array([0.1, 0.1, 0.1, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])


def integrated_score(S, W):
    """Compute integrated score with missing-data renormalization."""
    mask = ~np.isnan(S)
    if not mask.any():
        return 0.0
    S_filled = np.where(np.isnan(S), 0.0, S)
    w_masked = W[mask]
    return np.sum(S_filled[mask] * w_masked) / w_masked.sum()


def main():
    print("=== Leave-One-STREAM-Out Analysis (stream ablation) ===")

    p = RUN_DIR / "rank.csv"
    if not p.exists():
        print("Error: run the pipeline first (rank.csv missing)")
        return 1

    df = pd.read_csv(p).drop_duplicates(subset="gene_id", keep="first")

    gene_col = "gene_id" if "gene_id" in df.columns else "gene_id_v6"
    name_col = "gene_name" if "gene_name" in df.columns else gene_col
    score_col = "integrated_score" if "integrated_score" in df.columns else "composite_score"

    stream_cols = [s for s in STREAMS if s in df.columns]
    if len(stream_cols) < 3:
        stream_cols = [c for c in df.columns if c in STREAMS]

    scores_matrix = df[stream_cols].values.astype(float)
    # weights looked up BY STREAM NAME (never prefix-truncated: a missing
    # middle stream would silently misalign every subsequent weight)
    weights = np.array([W_DEFAULT[STREAMS.index(s)] for s in stream_cols])
    weights = weights / weights.sum()

    baseline_scores = np.array([integrated_score(row, weights) for row in scores_matrix])
    # 2026-09-06 audit fix (deterministic tie-broken ordering): the score
    # family has MASSIVE ties (541 genes at 0.5667, 425 at 0.6000, 1,817
    # in tie-groups >=5; four genes tied at 0.9667 spanning baseline
    # ranks 6-9). The previous np.argsort(-score) resolved top-10
    # boundary ties by ROW ORDER (rank.csv happens to be sorted by
    # descending integrated_score, so this silently favored genes
    # already higher in the baseline — biasing overlap toward stability
    # — and pd.Series.rank(method='min') produced corrupted rank
    # matrices with impossible repeated ranks and gaps). Both are now
    # tie-broken by gene_id ascending, matching select_top's discipline.
    gene_ids = df[gene_col].astype(str).values

    def _order_desc(scores: np.ndarray) -> np.ndarray:
        """Indices sorting scores descending, ties broken by gene_id
        ascending (lexicographic inversion of the id, same trick as
        prioritize.select_top)."""
        gid_desc = np.array(["".join(chr(0x10FFFF - ord(ch)) for ch in g)
                             for g in gene_ids])
        return np.lexsort((gid_desc, -scores))

    def _rank_series(scores: np.ndarray) -> np.ndarray:
        """Dense 1-based ranks with the same deterministic tie-break —
        identical genes get identical ranks, no gaps, boundary
        membership matches _order_desc exactly."""
        order = _order_desc(scores)
        ranks = np.empty(len(scores), dtype=int)
        # dense ranking over distinct score values; ties share the rank
        sorted_scores = scores[order]
        dense = np.empty(len(scores), dtype=int)
        dense[0] = 1
        for k in range(1, len(scores)):
            dense[k] = dense[k - 1] + (0 if sorted_scores[k] == sorted_scores[k - 1] else 1)
        ranks[order] = dense
        return ranks

    baseline_order = _order_desc(baseline_scores)
    top10_indices = baseline_order[:10]
    top10_ids = df[gene_col].values[top10_indices]
    top10_names = df[name_col].fillna(df[gene_col]).values[top10_indices]

    print(f"Candidates: {len(df)}, Streams: {len(stream_cols)}")
    print(f"Baseline top-10: {list(top10_names)}")

    results = []
    # Rank matrix for heatmap: rows = top10 candidates, cols = excluded streams
    rank_matrix_dict = {
        "gene_id": list(top10_ids),
        "gene_name": list(top10_names),
        "full_rank": list(range(1, 11)),
        "full_score": [float(baseline_scores[idx]) for idx in top10_indices],
    }

    for i, stream in enumerate(stream_cols):
        loo_streams = [s for j, s in enumerate(stream_cols) if j != i]
        loo_idx = [j for j in range(len(stream_cols)) if j != i]
        loo_weights = np.delete(weights, i)
        loo_weights = loo_weights / loo_weights.sum()
        loo_matrix = scores_matrix[:, loo_idx]

        loo_scores = np.array([integrated_score(row, loo_weights) for row in loo_matrix])
        loo_ranks = _rank_series(loo_scores)
        loo_order = _order_desc(loo_scores)
        loo_top10 = set(df[gene_col].values[loo_order[:10]])

        overlap = set(top10_ids) & loo_top10
        jaccard = len(overlap) / len(set(top10_ids) | loo_top10) if len(set(top10_ids) | loo_top10) > 0 else 0
        rank_corr = pd.Series(baseline_scores).corr(pd.Series(loo_scores), method="spearman")

        # Record ranks of baseline top-10 when this stream is excluded
        rank_matrix_dict[stream] = [int(loo_ranks[idx]) for idx in top10_indices]

        results.append({
            "excluded_stream": stream,
            "n_remaining_streams": len(loo_streams),
            "top10_overlap": len(overlap),
            "top10_jaccard": float(jaccard),
            "spearman_correlation": float(rank_corr),
            "overlap_genes": sorted(list(overlap)),
            "new_in_top10": sorted(list(loo_top10 - set(top10_ids))),
            "dropped_from_top10": sorted(list(set(top10_ids) - loo_top10)),
        })
        print(f"  Exclude {stream:>20s}: overlap={len(overlap)}/10, "
              f"Jaccard={jaccard:.3f}, rho={rank_corr:.4f}")

    # Write rank matrix to loo_atlas_stability.csv for figure 27
    # (filename kept for compatibility; content is stream-ablation)
    matrix_df = pd.DataFrame(rank_matrix_dict)
    out_path = RESULTS_DIR / "loo_atlas_stability.csv"
    matrix_df.to_csv(out_path, index=False)
    print(f"\nSaved stability matrix: {out_path}")

    # Write summary metrics to loo_atlas_summary.csv
    summary_df = pd.DataFrame(results)
    summary_path = RESULTS_DIR / "loo_atlas_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"Saved summary metrics: {summary_path}")

    avg_overlap = np.mean([r["top10_overlap"] for r in results])
    avg_jaccard = np.mean([r["top10_jaccard"] for r in results])
    avg_rho = np.mean([r["spearman_correlation"] for r in results])
    print(f"\nMean overlap: {avg_overlap:.1f}/10")
    print(f"Mean Jaccard: {avg_jaccard:.3f}")
    print(f"Mean Spearman rho: {avg_rho:.4f}")
    # 2026-09-06 honesty note: mean Spearman rho over 11,675 mostly-tied
    # bulk scores is INSENSITIVE to top-10 churn (observed rho>=0.98 for
    # every ablation while top-10 overlap swings 5-10/10). The top-10
    # overlap/Jaccard columns are the stability headline; rho is a bulk
    # descriptor only.
    print("(stability headline = top10_overlap/jaccard; mean rho is a "
          "bulk-score descriptor, insensitive to top-rank churn)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
