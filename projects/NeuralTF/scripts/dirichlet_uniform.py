#!/usr/bin/env python
"""Dirichlet-robust prioritization with NON-INFORMATIVE (uniform) prior.

Unlike `dirichlet_prioritize.py` (which centers the Dirichlet on the default
weights), this script uses a **uniform Dirichlet** — no prior preference for
any weighting. This answers: "If we knew nothing about the relative importance
of the 9 evidence streams, which candidates emerge as best supported across
ALL plausible weightings?"

Method:
    alpha_i = 1 for all 9 streams  (uniform over the 9-simplex)
    Sample 1000 weight vectors, apply each to ALL candidates,
    take the median integrated score.

SCIENTIFIC RATIONALE:
The centered Dirichlet (k=40 around defaults) is a sensitivity analysis:
"Assuming defaults are approximately right, how robust is the ranking?"
The uniform Dirichlet is a prior-free analysis:
"What does the data itself say when no weighting preference is imposed?"
A candidate that ranks high under BOTH is fundamentally robust.

Outputs (CSVs/TXT into `projects/NeuralTF/results/`, gitignored):
  - dirichlet_uniform_top10.csv            (track-based 5A+5B)
  - dirichlet_uniform_overall_top10.csv      (overall top-10 by uniform median)
  - dirichlet_uniform_full_rank.csv          (all 99 candidates with both scores)
  - dirichlet_uniform_summary.txt            (3-way comparison stats)

Usage:
    python projects/NeuralTF/scripts/dirichlet_uniform.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO = Path(__file__).resolve().parents[3]
DATA = REPO / "projects" / "NeuralTF" / "data"
RUN = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
OUT = REPO / "projects" / "NeuralTF" / "results"

# Import shared functions from the centered-Dirichlet script
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from dirichlet_prioritize import (
    STREAMS,
    N_DRAWS,
    SEED,
    integrated_scores,
    build_csv,
    clean_ortholog,
    rnai_marker_notes,
    read_mmc5,
)
from bioforge.projects.neuraltf.prioritize import attach_v4

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
# Uniform (non-informative) prior — all streams equal
ALPHA_UNIFORM = np.ones(len(STREAMS), dtype=float)

# Default weights (expression=0.2, all others=0.1)
W_DEFAULT = np.array([0.200, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100])


# ---------------------------------------------------------------------------
# Uniform Dirichlet sampling
# ---------------------------------------------------------------------------
def dirichlet_uniform_median_scores(S: np.ndarray,
                                    n_draws: int,
                                    seed: int) -> np.ndarray:
    """Median integrated score across draws under uniform Dirichlet prior.

    Each draw samples ONE weight vector uniformly from the 9-simplex
    (alpha_i = 1 for all streams). The SAME weight vector is applied to
    ALL candidates per draw.
    """
    rng = np.random.default_rng(seed)
    n = S.shape[0]
    all_scores = np.empty((n_draws, n), dtype=np.float32)
    mask = ~np.isnan(S)

    for d in range(n_draws):
        # Sample ONE weight vector uniformly from the simplex
        w = rng.dirichlet(ALPHA_UNIFORM)
        # Apply to all candidates (NaN streams zeroed)
        all_scores[d] = integrated_scores(S, w)

    return np.median(all_scores, axis=0)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print("== Dirichlet-UNIFORM (non-informative) prioritization ==")
    print(f"  Prior: Dirichlet(alpha_i=1) — uniform over 9-simplex")
    print(f"  Draws: {N_DRAWS}, seed: {SEED}")

    # Load data
    rank_path = RUN / "rank_neural.csv"
    rank = pd.read_csv(rank_path)
    print(f"  candidates: {len(rank)}")

    # Load bridge for v6→v4 mapping
    bridge_path = DATA / "bridge.csv"
    bridge = pd.read_csv(bridge_path, dtype=str)

    # Extract stream matrix
    S = rank[STREAMS].to_numpy(dtype=float)

    # Compute uniform-Dirichlet median scores
    median_scores = dirichlet_uniform_median_scores(S, N_DRAWS, SEED)
    rank["uniform_median_score"] = median_scores

    # Also compute fixed-weight scores for comparison
    fixed_mask = ~np.isnan(S)
    rank["fixed_score"] = np.where(fixed_mask, S, 0.0) @ W_DEFAULT

    print(f"  Uniform median scores: "
          f"min={median_scores.min():.4f}, max={median_scores.max():.4f}, "
          f"mean={median_scores.mean():.4f}")

    # --- 3-way comparison ---------------------------------------------------
    print("\n  === Overall top-10 by uniform-Dirichlet median ===")
    top_uniform = rank.nlargest(10, "uniform_median_score")
    for _, r in top_uniform.iterrows():
        print(f"    {r['gene_name']:>8}  uniform={r['uniform_median_score']:.4f}  "
              f"fixed={r['fixed_score']:.4f}  proof={r['proof_status']}")

    # Fixed-weight top-10 for comparison
    print("\n  === Fixed-weight overall top-10 ===")
    top_fixed = rank.nlargest(10, "fixed_score")
    for _, r in top_fixed.iterrows():
        print(f"    {r['gene_name']:>8}  fixed={r['fixed_score']:.4f}  "
              f"uniform={r['uniform_median_score']:.4f}")

    # Compare
    unif_set = set(top_uniform["gene_id"].tolist())
    fix_set = set(top_fixed["gene_id"].tolist())
    print(f"\n  Top-10 overlap (uniform vs fixed): {len(unif_set & fix_set)}/10")
    if unif_set != fix_set:
        print(f"  Fixed-only: {[rank[rank['gene_id']==g]['gene_name'].iloc[0] for g in fix_set - unif_set]}")
        print(f"  Uniform-only: {[rank[rank['gene_id']==g]['gene_name'].iloc[0] for g in unif_set - fix_set]}")
    else:
        print("  Top-10 identical under both methods")

    # --- Track-based top-10 (5A + 5B) --------------------------------------
    print("\n  === Track-based top-5 per track (uniform prior) ===")

    # Need annotations for Track B filter — load from PlanMine parquet
    ann_path = REPO / "datasets" / "processed" / "planmine_annotations.parquet"
    ann = pd.read_parquet(ann_path)
    # Collapse to one row per gene
    from bioforge.projects.neuraltf.prioritize import (
        summarize_annotations, prepare_candidates, assign_tracks,
        select_top, merge_annotations, map_v6_to_v4,
        compute_composite,
    )

    KING_DIR = REPO / "datasets" / "raw" / "Supplementary_Data_ King_2024"

    def _resolve(king_dir: Path, name: str) -> Path:
        cand = king_dir / f"1-s2.0-S2211124724001712-{name}.xlsx"
        if cand.exists():
            return cand
        if king_dir.exists():
            for p in sorted(king_dir.iterdir()):
                if p.suffix.lower() == ".xlsx" and p.stem.lower().endswith(name):
                    return p
        return cand

    def read_mmc4(path: Path) -> pd.DataFrame:
        raw = pd.read_excel(path, header=None)
        header_row = None
        for i in range(min(len(raw), 6)):
            vals = [str(x) for x in raw.iloc[i].tolist()[:8]]
            if "Gene ID" in vals and "Human Best Blast Hit" in vals:
                header_row = i
                break
        if header_row is None:
            return pd.DataFrame()
        df = pd.DataFrame(raw.iloc[header_row + 1:].values,
                          columns=raw.iloc[header_row].tolist())
        df = df.dropna(subset=["Gene ID"]).reset_index(drop=True)
        df["Gene ID"] = df["Gene ID"].astype(str).str.strip()
        return df

    mmc4 = None
    mmc4_path = _resolve(KING_DIR, "mmc4")
    if mmc4_path.exists():
        mmc4 = read_mmc4(mmc4_path)

    ann_sum = summarize_annotations(ann)
    cand = prepare_candidates(rank, mmc4=mmc4)
    cand["integrated_score"] = median_scores  # Use uniform median as base
    cand = merge_annotations(cand, ann_sum)
    cand = compute_composite(cand)

    a, b = assign_tracks(cand)
    b_filtered = b[
        (b["dna_binding_domains"].astype(str).str.strip() != "")
        | (b["mmc4_tf_flag"].astype(str).str.upper() == "TF")
    ]
    print(f"  Track B after TF-domain filter: {len(b_filtered)}/{len(b)}")

    ta = select_top(a, 5).assign(track="A")
    tb = select_top(b_filtered, 5).assign(track="B")
    top_track = pd.concat([ta, tb], ignore_index=True)

    # Attach v4 ID (matches dirichlet_prioritize.py build_csv expectations)
    mapping = map_v6_to_v4(bridge)
    top_track = attach_v4(top_track, mapping)

    # Rename for CSV consistency (matches dirichlet_prioritize.py pattern)
    top_track = top_track.rename(columns={"integrated_score": "dirichlet_median_score"})

    # --- Load mmc5 for RNAi notes (matching centered-Dirichlet output) -----
    mmc5 = None
    mmc5_path = _resolve(KING_DIR, "mmc5")
    if mmc5_path.exists():
        mmc5 = read_mmc5(mmc5_path)

    # RNAi phenotype notes
    notes = []
    for _, r in top_track.iterrows():
        if r["proof_status"] == "known_rnai_validated":
            notes.append(rnai_marker_notes(mmc5, r["gene_id"]))
        else:
            notes.append("Not RNAi-tested in King 2024 mmc5; novel neural-fate candidate")
    top_track["rnai_phenotype_notes"] = notes

    print("\n  Track A top-5:")
    for _, r in top_track[top_track["track"] == "A"].iterrows():
        print(f"    {r['gene_name']:>8}  composite={r['composite_score']:.3f}  "
              f"domain={(r.get('dna_binding_domains') or 'none')[:30]}")
    print("\n  Track B top-5:")
    for _, r in top_track[top_track["track"] == "B"].iterrows():
        print(f"    {r['gene_name']:>8}  composite={r['composite_score']:.3f}  "
              f"domain={(r.get('dna_binding_domains') or 'none')[:30]}")

    # --- Save CSVs ----------------------------------------------------------
    OUT.mkdir(parents=True, exist_ok=True)
    # Primary track-based CSV (same column structure as dirichlet_top10_prioritized.csv)
    csv_out = build_csv(top_track)
    # Rename dirichlet_median_score → uniform_median_score for clarity
    csv_out = csv_out.rename(columns={"dirichlet_median_score": "uniform_median_score"})
    csv_path = OUT / "dirichlet_uniform_top10.csv"
    csv_out.to_csv(csv_path, index=False)
    print(f"\n  wrote {csv_path.name} (track-based 5A+5B, {len(csv_out)} rows)")

    # Overall top-10 by uniform median
    overall_csv = pd.DataFrame({
        "gene_id_v6": top_uniform["gene_id"],
        "gene_name": top_uniform["gene_name"],
        "uniform_median_score": top_uniform["uniform_median_score"],
        "fixed_weight_score": top_uniform["fixed_score"],
        "proof_status": top_uniform["proof_status"],
    })
    overall_csv.to_csv(OUT / "dirichlet_uniform_overall_top10.csv", index=False)
    print(f"  wrote dirichlet_uniform_overall_top10.csv (overall top-10)")

    # Full-rank CSV: all 99 candidates with both uniform and fixed scores
    # (matches weight_sensitivity_draws.csv pattern — all candidates, all scores)
    full_rank = pd.DataFrame({
        "gene_id_v6": rank["gene_id"],
        "gene_name": rank["gene_name"],
        "uniform_median_score": rank["uniform_median_score"],
        "fixed_weight_score": rank["fixed_score"],
        "proof_status": rank["proof_status"],
    })
    full_rank_csv = OUT / "dirichlet_uniform_full_rank.csv"
    full_rank.to_csv(full_rank_csv, index=False)
    print(f"  wrote {full_rank_csv.name} (all {len(full_rank)} candidates)")

    # --- Load centered-Dirichlet results for 3-way comparison ---------------
    centered_csv = OUT / "dirichlet_top10_prioritized.csv"
    if centered_csv.exists():
        centered = pd.read_csv(centered_csv)
        print(f"\n  === 3-way comparison (fixed / centered / uniform) ===")
        # Merge all three scores
        comp = top_track.merge(
            centered[["gene_id_v6", "dirichlet_median_score"]].rename(
                columns={"dirichlet_median_score": "centered_median"}),
            left_on="gene_id", right_on="gene_id_v6", how="left"
        ).merge(
            top_track[["gene_id", "dirichlet_median_score"]].rename(
                columns={"dirichlet_median_score": "uniform_median"}),
            left_on="gene_id", right_on="gene_id", how="left"
        )

        # Compute scores from raw S for each
        for _, r in comp.iterrows():
            gid = r["gene_id"]
            row = rank[rank["gene_id"] == gid].iloc[0]
            fixed = row["fixed_score"]
            unif = row["uniform_median_score"]
            cent = r["centered_median"] if pd.notna(r["centered_median"]) else 0.0

            print(f"  {r['gene_name']:>8} (Track {r['track']}): "
                  f"fixed={fixed:.3f}  centered={cent:.3f}  uniform={unif:.3f}  "
                  f"range={max(fixed, cent, unif)-min(fixed, cent, unif):.3f}")

        # Summary stats
        fixed_vals = comp["gene_id"].apply(
            lambda g: rank[rank["gene_id"] == g].iloc[0]["fixed_score"])
        unif_vals = comp["uniform_median"]
        cent_vals = comp["centered_median"]

        max_range = max(
            abs(fixed_vals.dropna().max() - fixed_vals.dropna().min()) if len(fixed_vals.dropna()) > 0 else 0,
            abs(unif_vals.dropna().max() - unif_vals.dropna().min()) if len(unif_vals.dropna()) > 0 else 0,
            abs(cent_vals.dropna().max() - cent_vals.dropna().min()) if len(cent_vals.dropna()) > 0 else 0,
        )

        summary_lines = [
            "Dirichlet-UNIFORM Prior: Summary Report",
            "=" * 50,
            "",
            f"Method: Uniform Dirichlet prior (alpha_i=1) over 9-simplex",
            f"Draws: {N_DRAWS}, Seed: {SEED}",
            "",
            "Interpretation: Non-informative prior — no preference for any",
            "weighting. Tests fundamental robustness independent of defaults.",
            "",
            f"Overall top-10 (by uniform median): {len(top_uniform)} candidates",
            f"Fixed-weight top-10 overlap: {len(unif_set & fix_set)}/10",
            "",
            "3-way score comparison (track-based 5A+5B):",
            "-" * 50,
        ]
        for _, r in comp.iterrows():
            gid = r["gene_id"]
            row = rank[rank["gene_id"] == gid].iloc[0]
            summary_lines.append(
                f"  {r['gene_name']:>8} (Track {r['track']}): "
                f"fixed={row['fixed_score']:.4f}  "
                f"centered={r['centered_median']:.4f}  "
                f"uniform={r['uniform_median']:.4f}"
            )

        summary_lines.extend([
            "",
            "Key insight:",
            "- Small score range (max - min) = robust across all 3 methods",
            "- Large score range = candidate depends on weight choice",
        ])

        summary_path = OUT / "dirichlet_uniform_summary.txt"
        summary_path.write_text("\n".join(summary_lines), encoding="utf-8")
        print(f"\n  wrote {summary_path.name}")

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
