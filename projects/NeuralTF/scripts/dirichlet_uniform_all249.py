#!/usr/bin/env python
"""Dirichlet-uniform prioritization on ALL 249 TF candidates (no neural filter).

Unlike `dirichlet_uniform.py` (which scores only the 99 neural-filtered
candidates), this script runs the non-informative Dirichlet sampler on the
**full 249-candidate atlas** — i.e., every TF with significant differential
expression (p ≤ 0.05) in at least one atlas, regardless of neural signal.

SCIENTIFIC RATIONALE:
The 99-neural analysis asks "are the neural TF picks robust?". This 249-wide
analysis asks "are the best TFs overall robust, and do the neural ones
emerge without any neural filtering?". If the same candidates rise to the
top under both scopes, the neural filter is doing real work (selecting
biologically meaningful candidates), not just cherry-picking.

Method:
    alpha_i = 1 for all 7 streams  (uniform over the 7-simplex)
    Sample 1000 weight vectors, apply each to ALL 249 candidates
    Take the median integrated score across draws.

Track A/B selection still applies (RNAi-validated → A; TF-domain → B) but
across a wider pool. Candidates lacking neural evidence may still surface if
they have strong expression + reproducibility signals.

Outputs (CSVs in `projects/NeuralTF/results/`, gitignored):
  - dirichlet_uniform_all249_top10.csv           (track-based 5A+5B)
  - dirichlet_uniform_all249_overall_top10.csv   (overall top-10 by score)
  - dirichlet_uniform_all249_full_rank.csv       (all 249 with scores)
  - dirichlet_uniform_all249_summary.txt         (3-way: 99-fixed / 99-uniform / 249-uniform)

Usage:
    python projects/NeuralTF/scripts/dirichlet_uniform_all249.py
"""
from __future__ import annotations

from pathlib import Path
import sys

# Import shared functions from the centered-Dirichlet script
sys.path.insert(0, str(Path(__file__).resolve().parent))
from dirichlet_prioritize import (
    STREAMS,
    N_DRAWS,
    SEED,
    integrated_scores,
    build_csv,
    clean_ortholog,
    rnai_marker_notes,
)
from dirichlet_uniform import (
    ALPHA_UNIFORM,
    dirichlet_uniform_median_scores,
)
from bioforge.projects.neuraltf.prioritize import (
    attach_v4,
    assign_tracks,
    compute_composite,
    merge_annotations,
    prepare_candidates,
    select_top,
    summarize_annotations,
)

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO   = Path(__file__).resolve().parents[3]
DATA   = REPO / "projects" / "NeuralTF" / "data"
RUN    = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
OUT    = REPO / "projects" / "NeuralTF" / "results"
OUT.mkdir(parents=True, exist_ok=True)

KING_DIR = REPO / "datasets" / "raw" / "Supplementary_Data_ King_2024"


# ---------------------------------------------------------------------------
# Helpers (reused from dirichlet_prioritize.py)
# ---------------------------------------------------------------------------
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


def read_mmc5(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, header=None)
    body_start = 4
    for i in range(min(len(raw), 6)):
        vals = [str(x) for x in raw.iloc[i].tolist()[:4]]
        if any("FSTF" in v for v in vals):
            body_start = i + 1
            break
    df = raw.iloc[body_start:].dropna(how="all").reset_index(drop=True)
    df.columns = ["fstf_rnai"] + [f"marker_{j}" for j in range(1, df.shape[1])]
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print("== Dirichlet-UNIFORM (all-249) prioritization ==")
    print(f"  Prior: Dirichlet(alpha_i=1) — uniform over 7-simplex")
    print(f"  Draws: {N_DRAWS}, seed: {SEED}")

    # --- Load data ---------------------------------------------------------
    rank_path = RUN / "rank.csv"
    rank = pd.read_csv(rank_path)
    print(f"  candidates: {len(rank)} (full atlas, no neural filter)")

    bridge_path = DATA / "bridge.csv"
    bridge = pd.read_csv(bridge_path, dtype=str)
    mmc4_path = _resolve(KING_DIR, "mmc4")
    mmc4 = read_mmc4(mmc4_path) if mmc4_path.exists() else pd.DataFrame()
    mmc5_path = _resolve(KING_DIR, "mmc5")
    mmc5 = read_mmc5(mmc5_path) if mmc5_path.exists() else pd.DataFrame()

    ann_path = REPO / "datasets" / "processed" / "planmine_annotations.parquet"
    ann = pd.read_parquet(ann_path) if ann_path.exists() else pd.DataFrame()

    # --- Compute uniform-Dirichlet median scores ----------------------------
    S = rank[STREAMS].to_numpy(dtype=float)
    median_scores = dirichlet_uniform_median_scores(S, N_DRAWS, SEED)
    rank["uniform_median_score"] = median_scores

    print(f"  Uniform median scores: "
          f"min={median_scores.min():.4f}, max={median_scores.max():.4f}, "
          f"mean={median_scores.mean():.4f}")

    # --- Build candidate frame with composite bonuses ----------------------
    from bioforge.projects.neuraltf.prioritize import map_v6_to_v4
    mapping = map_v6_to_v4(bridge)
    ann_sum = summarize_annotations(ann) if not ann.empty else pd.DataFrame()
    cand = prepare_candidates(rank, mmc4=mmc4)
    cand = attach_v4(cand, mapping)
    if not ann_sum.empty:
        cand = merge_annotations(cand, ann_sum)
    cand["integrated_score"] = median_scores
    cand = compute_composite(cand)
    print(f"  v6->v4 mapping: {cand['v4_mapping_flag'].value_counts().to_dict()}")

    # --- Track assignment + top 5 per track --------------------------------
    a, b = assign_tracks(cand)
    b_filtered = b[
        (b["dna_binding_domains"].astype(str).str.strip() != "")
        | (b["mmc4_tf_flag"].astype(str).str.upper() == "TF")
    ]
    print(f"  Track B after TF-domain filter: {len(b_filtered)}/{len(b)}")

    ta = select_top(a, 5).assign(track="A")
    tb = select_top(b_filtered, 5).assign(track="B")
    top_track = pd.concat([ta, tb], ignore_index=True)
    top_track = top_track.rename(columns={"integrated_score": "dirichlet_median_score"})
    top_track = attach_v4(top_track, mapping)

    # RNAi phenotype notes
    notes = []
    for _, r in top_track.iterrows():
        if r["proof_status"] == "known_rnai_validated":
            notes.append(rnai_marker_notes(mmc5, r["gene_id"]))
        else:
            notes.append("Not RNAi-tested in King 2024 mmc5; novel candidate")
    top_track["rnai_phenotype_notes"] = notes

    print("\n  === Track A top-5 (all 249) ===")
    for _, r in top_track[top_track["track"] == "A"].iterrows():
        print(f"    {r['gene_name']:>8}  composite={r['composite_score']:.3f}  "
              f"domain={(r.get('dna_binding_domains') or 'none')[:30]}  "
              f"proof={r['proof_status']}")
    print("\n  === Track B top-5 (all 249) ===")
    for _, r in top_track[top_track["track"] == "B"].iterrows():
        print(f"    {r['gene_name']:>8}  composite={r['composite_score']:.3f}  "
              f"domain={(r.get('dna_binding_domains') or 'none')[:30]}  "
              f"proof={r['proof_status']}")

    # --- Save track-based CSV ----------------------------------------------
    csv_out = build_csv(top_track)
    csv_out = csv_out.rename(columns={"dirichlet_median_score": "uniform_median_score"})
    csv_path = OUT / "dirichlet_uniform_all249_top10.csv"
    csv_out.to_csv(csv_path, index=False)
    print(f"\n  wrote {csv_path.name} (track-based 5A+5B, {len(csv_out)} rows)")

    # --- Overall top-10 by uniform median ----------------------------------
    overall = rank.nlargest(10, "uniform_median_score").copy()
    overall_csv = pd.DataFrame({
        "gene_id_v6": overall["gene_id"],
        "gene_name": overall["gene_name"],
        "uniform_median_score": overall["uniform_median_score"],
        "proof_status": overall["proof_status"],
        "has_neural_signal": overall["neural_enriched"].notna(),
    })
    overall_csv.to_csv(OUT / "dirichlet_uniform_all249_overall_top10.csv", index=False)
    print(f"  wrote dirichlet_uniform_all249_overall_top10.csv "
          f"({len(overall_csv)} rows)")

    # --- Full-rank CSV (all 249) -------------------------------------------
    full_rank = pd.DataFrame({
        "gene_id_v6": rank["gene_id"],
        "gene_name": rank["gene_name"],
        "uniform_median_score": rank["uniform_median_score"],
        "proof_status": rank["proof_status"],
        "has_neural_signal": rank["neural_enriched"].notna(),
    })
    full_rank.to_csv(OUT / "dirichlet_uniform_all249_full_rank.csv", index=False)
    print(f"  wrote dirichlet_uniform_all249_full_rank.csv "
          f"(all {len(full_rank)} candidates)")

    # --- Summary comparison vs 99-neural uniform ---------------------------
    unif_99 = pd.read_csv(OUT / "dirichlet_uniform_top10.csv")
    unif_99_ids = set(unif_99["gene_id_v6"])
    all249_ids = set(csv_out["gene_id_v6"])
    overlap = unif_99_ids & all249_ids
    only_in_99 = unif_99_ids - all249_ids
    only_in_249 = all249_ids - unif_99_ids

    # How many top-10 candidates have neural signal?
    all249_with_neural = csv_out["gene_id_v6"].isin(
        rank[rank["neural_enriched"].notna()]["gene_id"]
    ).sum()

    summary = [
        "Dirichlet-UNIFORM (all-249) vs (99-neural) summary",
        "=" * 50,
        "",
        "Method: Uniform Dirichlet (alpha=1) prior.",
        f"  249-wide top-10 candidates: {len(all249_ids)}",
        f"  99-neural top-10 candidates: {len(unif_99_ids)}",
        f"  Overlap: {len(overlap)}/10",
        f"  Only in 99-neural: {sorted(only_in_99) if only_in_99 else 'none'}",
        f"  Only in 249-wide: {sorted(only_in_249) if only_in_249 else 'none'}",
        "",
        f"249-wide top-10 candidates with King neural signal: "
        f"{all249_with_neural}/10",
        "",
        "Interpretation:",
        "- High overlap (>=8/10) = neural filter is selecting robustly; "
        "the same candidates emerge regardless of scope",
        "- Low overlap (<=6/10) = the 99-neural filter is doing real work; "
        "different candidates win when the broader pool is included",
        f"- Current overlap = {len(overlap)}/10",
    ]
    summary_path = OUT / "dirichlet_uniform_all249_summary.txt"
    summary_path.write_text("\n".join(summary), encoding="utf-8")
    print(f"  wrote {summary_path.name}")
    print(f"\n  {len(overlap)}/10 overlap between 99-neural and 249-wide top-10")
    if only_in_99:
        print(f"  Only in 99-neural: {only_in_99}")
    if only_in_249:
        print(f"  Only in 249-wide: {only_in_249}")
    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
