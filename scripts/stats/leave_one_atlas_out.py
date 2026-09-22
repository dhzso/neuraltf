#!/usr/bin/env python
"""Leave-one-STREAM-out sensitivity analysis (dual-track shortlist + score-only views).

2026-09-04 audit relabel: despite its historical filename, this analysis
leaves out one EVIDENCE STREAM (of the 11), not one atlas - the atlas
contributions are already collapsed into single columns (expression is
max over Fincher/Plass/Cui/King), so individual-atlas attribution cannot
be isolated from rank.csv. The filename is kept for driver/figure
compatibility; the docstring, output columns, and figure labels say
"stream".

2026-09-19 rework (two metric families, tie-break bug fixed):

1. ``shortlist_*`` columns (the headline): each ablation RE-RUNS the
   published dual-track 5+5 selection (top-5 ``tested`` + top-5 gated
   ``not_tested`` by composite score, identical bonus mask and Track-B
   gate via ``bioforge.projects.neuraltf.prioritize``) on the ablated
   scores and reports overlap/turnover against the published shortlist
   (``top10_neural_tfs_prioritized.csv``). Previous runs measured a
   different set (the integrated-score top-10) and presented it as
   shortlist stability.

2. ``top10_*`` columns (legacy, kept for transparency): stability of the
   score-only top-10 (all candidates ranked by integrated score alone -
   no bonuses, no tracks, no gate). The historical ``_order_desc``
   sorted the inverted gene-id key ASCENDING, i.e. broke ties by gene_id
   DESCENDING - the opposite of ``prioritize.select_top``'s documented
   discipline. With 541-gene tie groups this could flip which genes sit
   at the top-10 boundary. Ties now break by gene_id ascending.

Usage:
    python scripts/stats/leave_one_atlas_out.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from bioforge.evidence.scoring import DEFAULT_WEIGHTS, STREAM_ORDER  # noqa: E402
from bioforge.projects.neuraltf.prioritize import (  # noqa: E402
    gate_track_b,
    select_top,
)

RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

STREAMS = [getattr(s, "value", s) for s in STREAM_ORDER]
W_BY_NAME = {getattr(k, "value", k): float(v)
             for k, v in DEFAULT_WEIGHTS.items()}


def integrated_score_rows(scores_matrix: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Vectorized per-record integrated score with missing-data renormalization.

    score(g) = sum(w_i * s_i over present streams) / sum(w_i over present),
    exactly mirroring bioforge.evidence.scoring.EvidenceScorer.
    """
    valid = ~np.isnan(scores_matrix)
    filled = np.nan_to_num(scores_matrix, nan=0.0)
    num = filled @ weights
    den = valid.astype(float) @ weights
    den = np.where(den > 0, den, 1.0)
    return num / den


def _select_shortlist(frame: pd.DataFrame) -> pd.DataFrame:
    """Run the published dual-track 5+5 selection on ``frame``.

    ``frame`` must carry: gene_id, integrated_score (the method base),
    n_streams, composite_score, proof_status, dna_binding_domains,
    mmc4_tf_flag. Uses the pipeline's own gate + select_top so the rule
    is identical to top10_neural_tfs_prioritized.csv.
    """
    df = frame.copy()
    df["composite_base_column"] = "integrated_score"
    a = df[df["proof_status"] == "tested"]
    b_mask = df["proof_status"] == "not_tested"
    b = df[b_mask]
    n_b_before = len(b)
    b = gate_track_b(b)
    if len(b) < 5:
        print(f"  [warn] Track-B pool after gate: {len(b)}/{n_b_before}")
    ta = select_top(a, 5).assign(track="A")
    tb = select_top(b, 5).assign(track="B")
    return pd.concat([ta, tb], ignore_index=True)


def main() -> int:
    print("=== Leave-One-STREAM-Out Analysis (stream ablation) ===")

    p = RUN_DIR / "rank.csv"
    if not p.exists():
        print("Error: run the pipeline first (rank.csv missing)")
        return 1

    # fixed_full_rank.csv is the deterministic full-candidate frame produced by
    # prioritize_neural_tfs.py: it carries the stream matrix AND the shared
    # bonus/gate/proof columns, so the ablated re-selection uses exactly the
    # published bonus mask and Track-B gate.
    ffr_path = RESULTS_DIR / "fixed_full_rank.csv"
    shortlist_path = RESULTS_DIR / "top10_neural_tfs_prioritized.csv"
    if not ffr_path.exists() or not shortlist_path.exists():
        print("Error: run scripts/prioritize_neural_tfs.py first "
              "(fixed_full_rank.csv / top10_neural_tfs_prioritized.csv missing)")
        return 1

    df = pd.read_csv(ffr_path).drop_duplicates(subset="gene_id", keep="first")
    published = pd.read_csv(shortlist_path)
    published = published.rename(columns={"gene_id_v6": "gene_id"})
    pub_ids = set(published["gene_id"])
    pub_by_track = {t: set(published.loc[published["track"] == t, "gene_id"])
                    for t in ("A", "B")}

    gene_col = "gene_id"
    name_col = "gene_name"
    stream_cols = [s for s in STREAMS if s in df.columns]
    if len(stream_cols) < 3:
        stream_cols = [c for c in df.columns if c in STREAMS]

    scores_matrix = df[stream_cols].to_numpy(dtype=float)
    # weights looked up BY STREAM NAME (never prefix-truncated: a missing
    # middle stream would silently misalign every subsequent weight)
    weights = np.array([W_BY_NAME[s] for s in stream_cols])
    weights = weights / weights.sum()

    gene_ids = df[gene_col].astype(str).values
    gene_names = df[name_col].fillna(df[gene_col]).astype(str).values
    bonus = pd.to_numeric(df["bonus_total"], errors="coerce").fillna(0.0).to_numpy()

    # ---- baseline (no ablation) -------------------------------------------
    baseline_scores = integrated_score_rows(scores_matrix, weights)
    n_streams_full = (~np.isnan(scores_matrix)).sum(axis=1)

    base_frame = pd.DataFrame({
        "gene_id": gene_ids,
        "integrated_score": baseline_scores,
        "n_streams": n_streams_full,
        "composite_score": baseline_scores + bonus,
        "proof_status": df["proof_status"].values,
        "dna_binding_domains": df["dna_binding_domains"].values,
        "mmc4_tf_flag": df["mmc4_tf_flag"].values,
    })
    baseline_shortlist = _select_shortlist(base_frame)
    base_ids = set(baseline_shortlist["gene_id"])
    if base_ids != pub_ids:
        print("  [warn] recomputed baseline shortlist != published shortlist:")
        print(f"    only-in-recomputed: {sorted(base_ids - pub_ids)}")
        print(f"    only-in-published : {sorted(pub_ids - base_ids)}")
    else:
        print("  baseline re-selection matches the published shortlist (10/10)")

    # ---- legacy score-only ordering (ties -> gene_id ASCENDING) -----------
    def _order_desc(scores: np.ndarray) -> np.ndarray:
        """Indices sorting scores descending, ties broken by gene_id
        ascending (np.lexsort sorts the last key first, both ascending:
        primary -scores, secondary gene_ids)."""
        return np.lexsort((gene_ids, -scores))

    def _rank_series(scores: np.ndarray) -> np.ndarray:
        """Dense 1-based ranks matching _order_desc's tie discipline."""
        order = _order_desc(scores)
        sorted_scores = scores[order]
        dense = np.empty(len(scores), dtype=int)
        dense[0] = 1
        for k in range(1, len(scores)):
            dense[k] = dense[k - 1] + (0 if sorted_scores[k] == sorted_scores[k - 1] else 1)
        ranks = np.empty(len(scores), dtype=int)
        ranks[order] = dense
        return ranks

    baseline_order = _order_desc(baseline_scores)
    score_top10_idx = baseline_order[:10]
    score_top10_ids = list(gene_ids[score_top10_idx])
    print(f"Candidates: {len(df)}, Streams: {len(stream_cols)}")
    print(f"Score-only top-10 (integrated score, no bonuses/tracks): "
          f"{[gene_names[i] for i in score_top10_idx]}")

    # ---- ablation loop -----------------------------------------------------
    results = []
    # stability matrix rows = the 10 PUBLISHED shortlist genes (dual-track view)
    pub_order = list(published["gene_id"])
    pub_rank_by_track = dict(zip(published["gene_id"], published["rank"]))
    pub_track = dict(zip(published["gene_id"], published["track"]))
    rank_matrix_dict = {
        "gene_id": pub_order,
        "gene_name": [str(published.loc[published["gene_id"] == g, "gene_name"].iloc[0])
                      for g in pub_order],
        "track": [pub_track[g] for g in pub_order],
        "full_rank": [pub_rank_by_track[g] for g in pub_order],
        "full_composite": [float(published.loc[published["gene_id"] == g,
                                               "composite_score"].iloc[0])
                           for g in pub_order],
    }
    idx_of = {g: i for i, g in enumerate(gene_ids)}

    for i, stream in enumerate(stream_cols):
        loo_idx = [j for j in range(len(stream_cols)) if j != i]
        loo_weights = np.delete(weights, i)
        loo_weights = loo_weights / loo_weights.sum()
        loo_matrix = scores_matrix[:, loo_idx]

        loo_scores = integrated_score_rows(loo_matrix, loo_weights)
        loo_n_streams = (~np.isnan(loo_matrix)).sum(axis=1)

        # (1) legacy score-only top-10 stability (fixed tie-break)
        loo_ranks = _rank_series(loo_scores)
        loo_order = _order_desc(loo_scores)
        loo_score_top10 = set(gene_ids[loo_order[:10]])
        overlap_score10 = set(score_top10_ids) & loo_score_top10
        union_score10 = set(score_top10_ids) | loo_score_top10
        rank_corr = pd.Series(baseline_scores).corr(pd.Series(loo_scores),
                                                    method="spearman")

        # (2) dual-track shortlist stability (the published 5+5 selection)
        abl_frame = pd.DataFrame({
            "gene_id": gene_ids,
            "integrated_score": loo_scores,
            "n_streams": loo_n_streams,
            "composite_score": loo_scores + bonus,
            "proof_status": df["proof_status"].values,
            "dna_binding_domains": df["dna_binding_domains"].values,
            "mmc4_tf_flag": df["mmc4_tf_flag"].values,
        })
        abl_shortlist = _select_shortlist(abl_frame)
        abl_ids = set(abl_shortlist["gene_id"])
        overlap_sl = pub_ids & abl_ids
        union_sl = pub_ids | abl_ids
        ta_overlap = len(pub_by_track["A"] & set(abl_shortlist.loc[abl_shortlist["track"] == "A", "gene_id"]))
        tb_overlap = len(pub_by_track["B"] & set(abl_shortlist.loc[abl_shortlist["track"] == "B", "gene_id"]))

        # stability matrix: published genes' within-track rank under ablation
        abl_rank = dict(zip(abl_shortlist["gene_id"], abl_shortlist["rank"]))
        rank_matrix_dict[stream] = [
            int(abl_rank[g]) if g in abl_rank else 0 for g in pub_order
        ]

        results.append({
            "excluded_stream": stream,
            "n_remaining_streams": len(stream_cols) - 1,
            "shortlist_overlap": len(overlap_sl),
            "shortlist_jaccard": len(overlap_sl) / len(union_sl) if union_sl else 0.0,
            "track_a_overlap": ta_overlap,
            "track_b_overlap": tb_overlap,
            "shortlist_new_in": sorted(list(abl_ids - pub_ids)),
            "shortlist_dropped": sorted(list(pub_ids - abl_ids)),
            "top10_overlap": len(overlap_score10),
            "top10_jaccard": len(overlap_score10) / len(union_score10) if union_score10 else 0.0,
            "spearman_correlation": float(rank_corr),
            "overlap_genes": sorted(list(overlap_score10)),
            "new_in_top10": sorted(list(loo_score_top10 - set(score_top10_ids))),
            "dropped_from_top10": sorted(list(set(score_top10_ids) - loo_score_top10)),
        })
        print(f"  Exclude {stream:>20s}: shortlist={len(overlap_sl)}/10 "
              f"(A {ta_overlap}/5, B {tb_overlap}/5), "
              f"score-top10={len(overlap_score10)}/10, rho={rank_corr:.4f}")

    # ---- outputs -----------------------------------------------------------
    matrix_df = pd.DataFrame(rank_matrix_dict)
    out_path = RESULTS_DIR / "loo_atlas_stability.csv"
    matrix_df.to_csv(out_path, index=False)
    print(f"\nSaved published-shortlist stability matrix: {out_path}")

    summary_df = pd.DataFrame(results)
    summary_path = RESULTS_DIR / "loo_atlas_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"Saved summary metrics: {summary_path}")

    sl_ov = [r["shortlist_overlap"] for r in results]
    t10_ov = [r["top10_overlap"] for r in results]
    avg_rho = np.mean([r["spearman_correlation"] for r in results])
    print(f"\nMean PUBLISHED-SHORTLIST overlap: {np.mean(sl_ov):.1f}/10  (headline)")
    print(f"Mean score-only top-10 overlap:   {np.mean(t10_ov):.1f}/10  (legacy view)")
    print(f"Mean Spearman rho: {avg_rho:.4f}")
    # honesty note: mean Spearman rho over 11,696 mostly-tied bulk scores is
    # INSENSITIVE to shortlist churn. The shortlist_overlap/jaccard columns
    # are the stability headline; rho is a bulk descriptor only.
    print("(stability headline = shortlist_overlap/jaccard; mean rho is a "
          "bulk-score descriptor, insensitive to top-rank churn)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
