#!/usr/bin/env python
"""Dirichlet-Centered robustness analysis on all candidates.

Evaluates ranking stability under Centered Dirichlet weight resampling
(k=40, 1,000 draws, seed=2024) across all candidates detected by the pipeline.

Unified method philosophy (WS2): identical candidate universe (rank.csv),
identical annotation summarization (one row per gene), identical bonus mask
and Track-B gate as the fixed method. The only difference between methods
is the weight vector used for the base score.

Outputs (in `projects/NeuralTF/results/`):
  - dirichlet_centered_full_rank.csv      — all candidates sorted by composite
  - dirichlet_centered_top10.csv          — 5 Track A + 5 Track B dual-track shortlist
  - dirichlet_centered_overall_top10.csv  — overall top-10 by composite score
  - dirichlet_centered_summary.txt        — summary of top ranked candidates
  - dirichlet_centered_draw_scores.csv    — per-candidate draw score matrix (CI inputs)

Usage:
    python projects/NeuralTF/scripts/dirichlet_centered.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Path setup
_script_dir = Path(__file__).resolve().parent
_root = _script_dir.parents[2]           # projects/NeuralTF/scripts -> repo root
sys.path.insert(0, str(_root / "src"))   # bioforge package

from bioforge.projects.neuraltf.prioritize import (  # noqa: E402
    apply_bonuses,
    assign_tracks,
    gate_track_b,
    map_v6_to_v4,
    merge_annotations,
    prepare_candidates,
    select_top,
    summarize_annotations,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
STREAMS = [
    "expression",
    "specificity",
    "reproducibility",
    "rnai",
    "correlation",
    "neural_enriched",
    "neural_specificity",
    "perez_lineage",
    "perez_influence",
]
W_DEFAULT = np.array([0.20, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10])
N_DRAWS = 1000
K_DIR = 40.0
SEED = 2024

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / "projects" / "NeuralTF" / "data"
RUN = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
OUT = REPO / "projects" / "NeuralTF" / "results"
OUT.mkdir(parents=True, exist_ok=True)


def _resolve(stem: str) -> Path | None:
    for base in (DATA, REPO / "datasets" / "processed", REPO / "datasets" / "raw"):
        for ext in (".parquet", ".csv", ".tsv"):
            p = base / f"{stem}{ext}"
            if p.exists():
                return p
    return None


def dirichlet_scores_all_draws(
    S: np.ndarray,
    W: np.ndarray,
    n_draws: int = N_DRAWS,
    k: float = K_DIR,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Draw Dirichlet weights centered at W; return the FULL (n_cand, n_draws)
    score matrix (per-candidate NaN-renormalized, mirroring the scorer)."""
    if rng is None:
        rng = np.random.default_rng(SEED)

    alpha = k * W
    draws = rng.dirichlet(alpha, size=n_draws)  # shape (n_draws, n_streams)

    valid_mask = ~np.isnan(S)                   # shape (n_cand, n_streams)
    S_filled = np.nan_to_num(S, nan=0.0)

    # scores_all: (n_cand, n_draws)
    num = S_filled @ draws.T                    # (n_cand, n_draws)
    den = valid_mask.astype(float) @ draws.T    # (n_cand, n_draws)
    den = np.where(den > 0, den, 1.0)
    return num / den


def build_csv(top: pd.DataFrame) -> pd.DataFrame:
    bridge_path = DATA / "bridge.csv"
    if bridge_path.exists():
        bridge = pd.read_csv(bridge_path, dtype=str)
        v6_col = "v6_id" if "v6_id" in bridge.columns else bridge.columns[1]
        v4_col = "v4_id" if "v4_id" in bridge.columns else bridge.columns[2]
        v4_map = dict(zip(bridge[v6_col].astype(str), bridge[v4_col].astype(str)))
    else:
        v4_map = {}

    out = top.copy()
    out["gene_id_v6"] = out["gene_id"]
    out["gene_id_v4"] = out["gene_id_v6"].map(v4_map).fillna("")

    gene_names = out.get("gene_name", out["gene_id_v6"])
    out["gene_symbol"] = gene_names.fillna(out["gene_id_v6"])

    def _notes(r: pd.Series) -> str:
        parts: list[str] = []
        if str(r.get("proof_status", "")).strip() == "known_rnai_validated":
            parts.append("Known validated neural TF (positive control)")
        orth = str(r.get("human_ortholog", "") or "").strip()
        if orth and orth.lower() not in ("nan", "none"):
            parts.append(f"Human TF ortholog: {orth}")
        if not parts:
            parts.append("Novel candidate identified by multi-atlas integration")
        return "; ".join(parts)

    out["rnai_screen_or_marker_notes"] = out.apply(_notes, axis=1)

    preferred_cols = [
        "track",
        "rank_within_track",
        "gene_id_v6",
        "gene_id_v4",
        "gene_symbol",
        "human_ortholog",
        "composite_score",
        "dirichlet_median_score",
        "bonus_total",
        "proof_status",
        "rnai_screen_or_marker_notes",
    ]
    cols = [c for c in preferred_cols if c in out.columns]
    return out[cols]


def main() -> int:
    print("=== Dirichlet Centered Robustness Analysis (All Candidates) ===")

    rank_path = RUN / "rank.csv"
    if not rank_path.exists():
        raise FileNotFoundError(f"Pipeline rank file not found: {rank_path}")

    raw_rank = pd.read_csv(rank_path)
    print(f"Loaded {len(raw_rank)} candidates from {rank_path.name}")

    # Annotations: collapse long -> one row per gene BEFORE merging.
    # (Raw long joins explode 278 rows to thousands of duplicated rows and
    # corrupt every downstream top-10 selection.)
    ann_path = _resolve("planmine_annotations")
    if ann_path is not None:
        annot = pd.read_parquet(ann_path) if ann_path.suffix == ".parquet" \
            else pd.read_csv(ann_path)
        annot = summarize_annotations(annot)
    else:
        annot = pd.DataFrame()

    cand = prepare_candidates(raw_rank)
    cand = merge_annotations(cand, annot)

    # Uniqueness guard: exactly one row per gene after every merge.
    assert cand["gene_id"].is_unique, (
        f"row explosion detected: {len(cand)} rows / "
        f"{cand['gene_id'].nunique()} unique genes"
    )

    # Streams
    available_streams = [s for s in STREAMS if s in cand.columns]
    w = W_DEFAULT[:len(available_streams)]
    w = w / w.sum()

    S = cand[available_streams].to_numpy(dtype=float)
    rng = np.random.default_rng(SEED)

    print(f"Running Centered Dirichlet (k={K_DIR}, n_draws={N_DRAWS}, "
          f"{len(available_streams)} streams)...")
    scores_all = dirichlet_scores_all_draws(S, w, N_DRAWS, K_DIR, rng)
    median_scores = np.median(scores_all, axis=1)
    lo = np.percentile(scores_all, 2.5, axis=1)
    hi = np.percentile(scores_all, 97.5, axis=1)
    cand["dirichlet_median_score"] = median_scores
    cand["dirichlet_ci95_lo"] = lo
    cand["dirichlet_ci95_hi"] = hi

    # Same bonus mask as the fixed method (WS2 unification).
    cand = apply_bonuses(cand, "dirichlet_median_score")
    cand["rank"] = cand["composite_score"].rank(ascending=False, method="min").astype(int)

    # Save full rank
    full_rank = cand.sort_values("composite_score", ascending=False).reset_index(drop=True)
    full_rank_path = OUT / "dirichlet_centered_full_rank.csv"
    full_rank.to_csv(full_rank_path, index=False)
    print(f"Saved full ranking ({len(full_rank)} candidates, "
          f"{full_rank['gene_id'].nunique()} unique genes): {full_rank_path}")

    # Per-draw score matrix (drives bootstrap CIs + convergence figures)
    draw_df = pd.DataFrame(
        np.column_stack([cand["gene_id"].to_numpy(), scores_all]),
        columns=["gene_id"] + [f"draw_{i}" for i in range(scores_all.shape[1])],
    )
    draw_path = OUT / "dirichlet_centered_draw_scores.csv"
    draw_df.to_csv(draw_path, index=False)
    print(f"Saved draw-score matrix ({scores_all.shape[0]} cand x "
          f"{scores_all.shape[1]} draws): {draw_path}")

    # Dual-track top 10 selection (5 Track A + 5 Track B) — same Track B
    # domain gate as the fixed method.
    track_a, track_b = assign_tracks(cand)
    track_b = gate_track_b(track_b)
    print(f"  Track B after TF-domain gate: {len(track_b)}")

    top_a = select_top(track_a, n=5)
    top_a["track"] = "A"
    top_a["rank_within_track"] = range(1, len(top_a) + 1)

    top_b = select_top(track_b, n=5)
    top_b["track"] = "B"
    top_b["rank_within_track"] = range(1, len(top_b) + 1)

    top10 = pd.concat([top_a, top_b], ignore_index=True)
    assert top10["gene_id"].nunique() == len(top10), (
        "top-10 shortlist contains duplicate genes"
    )
    overall_top10 = full_rank.head(10)

    top10_csv = build_csv(top10)
    top10_path = OUT / "dirichlet_centered_top10.csv"
    top10_csv.to_csv(top10_path, index=False)
    print(f"Saved top-10 shortlist (5 Track A + 5 Track B): {top10_path}")

    overall_path = OUT / "dirichlet_centered_overall_top10.csv"
    overall_csv = build_csv(overall_top10)
    overall_csv.to_csv(overall_path, index=False)
    print(f"Saved overall top-10: {overall_path}")

    # Summary text
    summary_path = OUT / "dirichlet_centered_summary.txt"
    lines = [
        "DIRICHLET CENTERED (k=40) ROBUSTNESS SUMMARY",
        f"Total candidates evaluated: {len(full_rank)}",
        f"Streams used: {len(available_streams)}",
        "",
        "Top-10 Shortlist (5 Track A + 5 Track B):",
    ]
    for _, r in top10_csv.iterrows():
        lines.append(f"  [{r.get('track','?')}] {r.get('gene_symbol','?'):<12} "
                     f"(v6: {r.get('gene_id_v6','?')}) "
                     f"composite={r.get('composite_score', 0):.4f} "
                     f"(median={r.get('dirichlet_median_score', 0):.4f}, "
                     f"bonus={r.get('bonus_total', 0):.2f})")
    lines.append("")
    lines.append("Overall Top-10 by Composite Score (median + bonuses):")
    for _, r in overall_csv.iterrows():
        lines.append(f"  {r.get('gene_symbol','?'):<12} (v6: {r.get('gene_id_v6','?')}) "
                     f"composite={r.get('composite_score', 0):.4f} "
                     f"(median={r.get('dirichlet_median_score', 0):.4f})")

    with open(summary_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved summary: {summary_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
