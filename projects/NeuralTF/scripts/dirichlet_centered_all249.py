#!/usr/bin/env python
"""Dirichlet-Centered robustness analysis on ALL TF candidates (n=249+).

Runs the same Centered Dirichlet (k=40, 1,000 draws, seed 2024) analysis
as `dirichlet_prioritize.py`, but on the **full candidate set** (`rank.csv`,
~249+ candidates) rather than only the 99 neural-filtered candidates.

SCIENTIFIC RATIONALE:
`dirichlet_prioritize.py` (centered, k=40) asks:
  "Are the top neural TF picks robust to small weight perturbations?"

This script adds the complementary question:
  "Across all atlas-detected candidates (neural + non-neural), do the same
   top TFs emerge under centered Dirichlet perturbation?"

If high-scoring candidates in the centered_all249 run match those in the
centered_neural99 run, it confirms the neural filter is enriching true
signal, not artificially inflating rank positions.

Outputs (in `projects/NeuralTF/results/`):
  dirichlet_centered_all249_full_rank.csv      — all candidates, sorted by Dirichlet median
  dirichlet_centered_all249_top10.csv          — 5 Track A + 5 Track B
  dirichlet_centered_all249_overall_top10.csv  — overall top-10 by Dirichlet median
  dirichlet_centered_all249_summary.txt        — 3-way comparison: centered99 / uniform249 / centered249

Usage:
    python projects/NeuralTF/scripts/dirichlet_centered_all249.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Reuse shared logic from the existing scripts
_script_dir = Path(__file__).resolve().parent
_root = _script_dir.parents[2]           # projects/NeuralTF/scripts -> repo root
sys.path.insert(0, str(_root / "src"))   # bioforge package
sys.path.insert(0, str(_script_dir))     # dirichlet_prioritize, dirichlet_uniform


from dirichlet_prioritize import (   # noqa: E402
    STREAMS,
    N_DRAWS,
    K_DIR,
    SEED,
    W_DEFAULT,
    _resolve,
    read_mmc4,
    read_mmc5,
    dirichlet_median_scores,
    build_csv,
    clean_ortholog,
    rnai_marker_notes,
)
from bioforge.projects.neuraltf.prioritize import (   # noqa: E402
    attach_v4,
    assign_tracks,
    compute_composite,
    merge_annotations,
    prepare_candidates,
    select_top,
    summarize_annotations,
    map_v6_to_v4,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO     = Path(__file__).resolve().parents[3]
DATA     = REPO / "projects" / "NeuralTF" / "data"
RUN      = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
OUT      = REPO / "projects" / "NeuralTF" / "results"
KING_DIR = REPO / "datasets" / "raw" / "Supplementary_Data_ King_2024"
OUT.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Helpers (mirrors dirichlet_uniform_all249.py)
# ---------------------------------------------------------------------------

def _rank_all_path() -> Path:
    """Prefer the renamed rank_all_candidates.csv; fall back to rank.csv."""
    for name in ("rank_all_candidates.csv", "rank.csv"):
        p = RUN / name
        if p.exists():
            return p
    raise FileNotFoundError(
        f"Neither rank_all_candidates.csv nor rank.csv found in {RUN}. "
        "Run the pipeline (scripts/run.py) first."
    )


def _build_report(
    top10: pd.DataFrame,
    overall_top10: pd.DataFrame,
    n_all: int,
    n_draws: int,
    k: float,
    median_scores: np.ndarray,
) -> str:
    lines = [
        "# Dirichlet-Centered (k=40) Robustness — All Candidates\n",
        f"Candidates evaluated: **{n_all}** (full rank.csv, no neural pre-filter).\n",
        f"Dirichlet draws: {n_draws}, concentration k={k}, seed {SEED}.\n",
        "## Overall Top-10 (by Dirichlet median score)\n",
        "| rank | gene_id | gene_name | dirichlet_median | composite |",
        "|------|---------|-----------|-----------------|-----------|",
    ]
    for i, (_, r) in enumerate(overall_top10.iterrows(), 1):
        lines.append(
            f"| {i} | {r['gene_id']} | {r.get('gene_name', '')} | "
            f"{r.get('dirichlet_median_score', 0):.4f} | "
            f"{r.get('composite_score', 0):.4f} |"
        )
    lines += [
        "",
        "## Score Statistics\n",
        f"- Min median score: {float(median_scores.min()):.4f}",
        f"- Max median score: {float(median_scores.max()):.4f}",
        f"- Mean median score: {float(median_scores.mean()):.4f}",
        f"- Std median score: {float(median_scores.std()):.4f}",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    rank_path = _rank_all_path()
    bridge_path = DATA / "bridge.csv"
    ann_path = REPO / "datasets" / "processed" / "planmine_annotations.parquet"
    mmc4_path = _resolve(KING_DIR, "mmc4")
    mmc5_path = _resolve(KING_DIR, "mmc5")

    print("== Dirichlet-Centered (k=40) — All Candidates ==")
    print(f"  rank       : {rank_path}")
    print(f"  annotations: {ann_path}")

    # ---------- Load inputs --------------------------------------------------
    rank = pd.read_csv(rank_path)
    bridge = pd.read_csv(bridge_path, dtype=str)
    mmc4 = read_mmc4(mmc4_path)
    mmc5 = read_mmc5(mmc5_path)
    ann = pd.read_parquet(ann_path)

    n_all = len(rank)
    print(f"  candidates : {n_all}")
    assert n_all > 0, "rank.csv is empty"

    # ---------- Build candidate frame ----------------------------------------
    mapping = map_v6_to_v4(bridge)
    ann_sum = summarize_annotations(ann)
    cand = prepare_candidates(rank, mmc4=mmc4)
    cand = attach_v4(cand, mapping)
    cand = merge_annotations(cand, ann_sum)
    print(f"  v6->v4 mapping: {cand['v4_mapping_flag'].value_counts().to_dict()}")

    # ---------- Dirichlet Centered (k=40) sampling ---------------------------
    S = cand[STREAMS].to_numpy(dtype=float)
    rng = np.random.default_rng(SEED)
    print(f"  Running Centered Dirichlet (k={K_DIR}, {N_DRAWS} draws) on {n_all} candidates...")
    median_scores = dirichlet_median_scores(S, W_DEFAULT, N_DRAWS, K_DIR, rng)

    cand["dirichlet_median_score"] = median_scores
    cand["integrated_score"] = median_scores  # use as ranking signal
    cand = compute_composite(cand)

    print(
        f"  Scores: min={median_scores.min():.4f}, "
        f"max={median_scores.max():.4f}, mean={median_scores.mean():.4f}"
    )

    # ---------- Full rank CSV ------------------------------------------------
    full_rank = cand.sort_values("dirichlet_median_score", ascending=False).reset_index(drop=True)
    full_rank["rank_centered_all"] = range(1, len(full_rank) + 1)
    full_rank_path = OUT / "dirichlet_centered_all249_full_rank.csv"
    full_rank.to_csv(full_rank_path, index=False)
    print(f"  Full rank written: {full_rank_path} ({len(full_rank)} rows)")

    assert len(full_rank) >= 99, (
        f"Expected >= 99 rows in full rank (neural set), got {len(full_rank)}"
    )

    # ---------- Overall top-10 -----------------------------------------------
    overall_top10 = full_rank.head(10)[
        ["gene_id", "gene_name", "dirichlet_median_score", "composite_score",
         "proof_status", "tier"]
    ].reset_index(drop=True)
    ot10_path = OUT / "dirichlet_centered_all249_overall_top10.csv"
    overall_top10.to_csv(ot10_path, index=False)
    print(f"  Overall top-10 written: {ot10_path}")

    # ---------- Track-based top-10 (5A + 5B) ---------------------------------
    a_track, b_track = assign_tracks(cand)
    b_track = b_track[
        (b_track["dna_binding_domains"].astype(str).str.strip() != "")
        | (b_track["mmc4_tf_flag"].astype(str).str.upper() == "TF")
    ]
    ta = select_top(a_track, 5).assign(track="A")
    tb = select_top(b_track, 5).assign(track="B")
    top10 = pd.concat([ta, tb], ignore_index=True)
    notes = []
    for _, r in top10.iterrows():
        if r["proof_status"] == "known_rnai_validated":
            notes.append(rnai_marker_notes(mmc5, r["gene_id"]))
        else:
            notes.append("Not RNAi-tested in King 2024 mmc5; novel candidate")
    top10["rnai_phenotype_notes"] = notes
    top10_csv = build_csv(top10)
    top10_path = OUT / "dirichlet_centered_all249_top10.csv"
    top10_csv.to_csv(top10_path, index=False)
    print(f"  Track top-10 (5A+5B) written: {top10_path}")

    # ---------- 3-way summary text -------------------------------------------
    neural_path = OUT / "dirichlet_centered_full_rank.csv"
    uniform249_path = OUT / "dirichlet_uniform_all249_full_rank.csv"
    lines = [
        f"== Dirichlet Robustness 3-Way Comparison ==",
        f"",
        f"Centered k=40, neural-filtered (rank_neural.csv):  {neural_path.name if neural_path.exists() else 'MISSING'}",
        f"Uniform  α=1,  all candidates (rank.csv):          {uniform249_path.name if uniform249_path.exists() else 'MISSING'}",
        f"Centered k=40, all candidates (rank.csv):          {full_rank_path.name}",
        f"",
        f"Full set ({n_all} candidates) top-10 by Centered Dirichlet median:",
    ]
    for i, (_, r) in enumerate(overall_top10.iterrows(), 1):
        lines.append(
            f"  {i:2d}. {(r.get('gene_name') or r['gene_id']):<20} "
            f"score={r['dirichlet_median_score']:.4f}  {r.get('proof_status','')}"
        )
    summary_path = OUT / "dirichlet_centered_all249_summary.txt"
    summary_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  Summary written: {summary_path}")

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
