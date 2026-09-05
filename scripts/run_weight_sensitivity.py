#!/usr/bin/env python
"""Weight-sensitivity analysis: 1000 random weight draws over ALL candidates.

2026-09-04 redesign:

- Universe: rank.csv (all candidates), matching the three published
  prioritization methods. The previous version used only the 128-gene
  neural subset, so its "top-10 stability" claims described a different
  universe than the published shortlist.
- Membership definition: the same composite dual-track scheme the
  published methods use (base + bonus mask, Track A/B split, Track-B
  gate, deterministic tie-breaks). The previous version compared a raw
  weighted-score top-10 (no bonuses, no gate) against the published
  composite shortlist - an apples-to-oranges baseline.
- Challenger attribution: when a baseline shortlist gene drops out of a
  draw's shortlist, the displacer recorded is the gene that actually took
  the vacated SLOT (track-aware), not the #1 gene of the draw.
- Determinism: boundary ties are resolved by the same tie-break column
  set used by select_top (score -> integrated -> n_streams -> gene_id
  ascending), never by unstable argsort artifact.

Draw weights come from a symmetric Dirichlet(1,...,1) (uniform over the
weight simplex); per-draw scores renormalize over each candidate's
available streams, mirroring EvidenceScorer exactly.

Outputs (written to projects/NeuralTF/figures/ for figures 06/07):
- weight_sensitivity_draws.csv: draw, gene_id, gene_name, rank, in_top_10
- weight_sensitivity_top10_challengers.csv: per-candidate summary
  (baseline_rank, frac_draws_in_top10, best/median rank, honest
  track-aware displacement info)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bioforge.projects.neuraltf.prioritize import (  # noqa: E402
    apply_bonuses,
    assign_tracks,
    gate_track_b,
    merge_annotations,
    prepare_candidates,
    summarize_annotations,
)

RUN_DIR = ROOT / "projects" / "NeuralTF" / "runs" / "pipeline_run"
FIG_DIR = ROOT / "projects" / "NeuralTF" / "figures"
RESULTS_DIR = ROOT / "projects" / "NeuralTF" / "results"
DATA_DIR = ROOT / "projects" / "NeuralTF" / "data"
RAW_DIR = ROOT / "datasets" / "raw"

STREAMS = ["expression", "specificity", "reproducibility", "rnai",
           "correlation", "neural_enriched", "neural_specificity",
           "perez_lineage", "perez_influence"]
N_DRAWS = 1000
SEED = 42
TIE_COLS = ["integrated_score", "n_streams"]


def weighted_scores(S: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Row-wise weighted scores with renormalization over available streams.

    Identical to EvidenceScorer.integrated_score: a missing stream never
    penalizes a candidate (absence of evidence is not evidence of absence).
    """
    mask = ~np.isnan(S)
    filled = np.where(mask, S, 0.0)
    denom = (mask * w).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        num = filled @ w
        return np.where(denom > 0, num / denom, 0.0)


def _gene_id_desc(gid) -> str:
    return "".join(chr(0x10FFFF - ord(ch)) for ch in str(gid))


def _load_annotations() -> pd.DataFrame:
    for base in (ROOT / "datasets" / "processed", DATA_DIR, RAW_DIR):
        for ext in (".parquet", ".csv", ".tsv"):
            p = base / f"planmine_annotations{ext}"
            if p.exists():
                long_df = pd.read_parquet(p) if ext == ".parquet" else pd.read_csv(p)
                return summarize_annotations(long_df)
    return pd.DataFrame()


def _load_mmc4() -> pd.DataFrame | None:
    king = RAW_DIR / "Supplementary_Data_ King_2024"
    if not king.exists():
        return None
    p = king / "1-s2.0-S2211124724001712-mmc4.xlsx"
    if p.exists():
        return pd.read_excel(p, sheet_name="TF")
    for q in sorted(king.iterdir()):
        if q.suffix.lower() == ".xlsx" and q.stem.lower().endswith("mmc4"):
            return pd.read_excel(q, sheet_name="TF")
    return None


def build_shortlist(cand: pd.DataFrame, scores: np.ndarray) -> pd.DataFrame:
    """One draw's composite dual-track 5+5 shortlist, mirroring the
    published methods: bonus mask on top of the draw's base score,
    Track A/B split, Track-B domain gate, deterministic tie-breaks."""
    d = cand.copy()
    d["_draw_base"] = scores
    d = apply_bonuses(d, "_draw_base")
    a, b = assign_tracks(d)
    b = gate_track_b(b)

    def top5(track: pd.DataFrame) -> pd.DataFrame:
        if track.empty:
            return track
        tmp = track.copy()
        tmp["_gid_desc"] = [_gene_id_desc(g) for g in tmp["gene_id"]]
        cols = ["composite_score"] + [c for c in TIE_COLS if c in tmp.columns] + ["_gid_desc"]
        return tmp.sort_values(cols, ascending=False).head(5)

    return pd.concat([top5(a), top5(b)], ignore_index=True)


def main() -> int:
    rank_path = RUN_DIR / "rank.csv"
    baseline_path = RESULTS_DIR / "top10_neural_tfs_prioritized.csv"
    if not rank_path.exists():
        print(f"Missing {rank_path} - run the pipeline first")
        return 1

    raw_rank = pd.read_csv(rank_path)
    annot = _load_annotations()
    cand = prepare_candidates(raw_rank, mmc4=_load_mmc4())
    cand = merge_annotations(cand, annot)
    assert cand["gene_id"].is_unique, "row explosion in sensitivity universe"

    available = [s for s in STREAMS if s in cand.columns]
    S = cand[available].to_numpy(dtype=float)
    n = len(cand)
    print(f"{n} candidates (full rank.csv universe), {N_DRAWS} weight draws")

    rng = np.random.default_rng(SEED)
    W = rng.dirichlet(np.ones(len(available)), size=N_DRAWS)  # uniform simplex

    # Baseline shortlist: the published fixed-method 5+5 (composite,
    # gated, tie-broken) - NOT a raw integrated_score top-10.
    if baseline_path.exists():
        baseline_df = pd.read_csv(baseline_path)
        gid_col = "gene_id_v6" if "gene_id_v6" in baseline_df.columns else "gene_id"
        baseline_shortlist = set(baseline_df[gid_col].astype(str))
    else:
        baseline_shortlist = build_shortlist(cand, weighted_scores(S, np.full(len(available), 1 / len(available))))
        baseline_shortlist = set(baseline_shortlist["gene_id"])
    print(f"Baseline shortlist: {len(baseline_shortlist)} genes")

    baseline_rank_map = dict(zip(
        cand["gene_id"],
        pd.to_numeric(cand["integrated_score"], errors="coerce").fillna(0)
        .rank(ascending=False, method="min").astype(int)))

    draw_rows = []
    in_top10_count = {g: 0 for g in cand["gene_id"]}
    best_rank = {g: n + 1 for g in cand["gene_id"]}
    rank_samples = {g: [] for g in cand["gene_id"]}
    displaced_by = {g: {} for g in cand["gene_id"]}

    for d in range(N_DRAWS):
        scores = weighted_scores(S, W[d])
        short = build_shortlist(cand, scores)

        # Full-universe rank under this draw (for rank distributions).
        r = pd.DataFrame({"gene_id": cand["gene_id"], "s": scores})
        r = r.assign(_gid_desc=[_gene_id_desc(g) for g in r["gene_id"]])
        order = r.sort_values(["s", "_gid_desc"], ascending=False).index.to_numpy()
        ranks = np.empty(n, dtype=int)
        ranks[order] = np.arange(1, n + 1)

        short_set = set(short["gene_id"])
        # Slot-aware challenger attribution: when a baseline gene misses the
        # draw's shortlist, the displacer is the draw-NEW gene that occupied
        # the same track's tail slot, not the draw's #1.
        new_in = short[~short["gene_id"].isin(baseline_shortlist)]
        dropped = baseline_shortlist - short_set
        entrants = short_set | set(new_in["gene_id"])

        for i, gid in enumerate(cand["gene_id"]):
            in10 = gid in short_set
            # Only persist rows for genes that ever enter a shortlist or
            # hold a top-30 rank - keeps the draws CSV at a readable size
            # instead of 1000 x 11,672 rows (420 MB) of mostly rank >1000.
            if in10 or gid in entrants or ranks[i] <= 30:
                draw_rows.append((d + 1, gid, cand["gene_name"].iloc[i], int(ranks[i]), in10))
            rank_samples[gid].append(int(ranks[i]))
            if int(ranks[i]) < best_rank[gid]:
                best_rank[gid] = int(ranks[i])
            if in10:
                in_top10_count[gid] += 1
            elif gid in dropped:
                # The displacer is the new gene in the same TRACK that
                # took the vacated slot (the last new entrant of that
                # track), not the draw's #1 gene.
                track_of = dict(zip(short["gene_id"], short.get("track", "?")))
                gtrack = track_of.get(gid, "A")
                if "track" in new_in.columns:
                    same_track_new = new_in[new_in["track"] == gtrack]
                else:
                    same_track_new = new_in
                displacer_pool = same_track_new if len(same_track_new) else new_in
                if len(displacer_pool):
                    displacer = str(displacer_pool["gene_id"].iloc[-1])
                    displaced_by[gid][displacer] = displaced_by[gid].get(displacer, 0) + 1

        if (d + 1) % 100 == 0:
            print(f"  {d + 1}/{N_DRAWS}", flush=True)

    draws_df = pd.DataFrame(draw_rows, columns=["draw", "gene_id", "gene_name", "rank", "in_top_10"])

    challenge_rows = []
    for i, gid in enumerate(cand["gene_id"]):
        if in_top10_count[gid] == 0:
            continue
        displacers = displaced_by.get(gid, {})
        most_displacer = max(displacers, key=displacers.get) if displacers else ""
        challenge_rows.append({
            "gene_id": gid,
            "gene_name": cand["gene_name"].iloc[i],
            "baseline_rank": int(baseline_rank_map.get(gid, n + 1)),
            "in_baseline_top10": gid in baseline_shortlist,
            "n_draws_in_top10": in_top10_count[gid],
            "frac_draws_in_top10": round(in_top10_count[gid] / N_DRAWS, 4),
            "best_rank": best_rank[gid],
            "median_rank": int(np.median(rank_samples[gid])),
            "displaced_most_often_by": most_displacer,
            "n_times_displaced": displacers.get(most_displacer, 0) if displacers else 0,
        })
    challengers_df = pd.DataFrame(challenge_rows).sort_values(
        "frac_draws_in_top10", ascending=False).reset_index(drop=True)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    draws_path = FIG_DIR / "weight_sensitivity_draws.csv"
    challengers_path = FIG_DIR / "weight_sensitivity_top10_challengers.csv"
    draws_df.to_csv(draws_path, index=False)
    challengers_df.to_csv(challengers_path, index=False)
    print(f"\nSaved {len(draws_df)} draw rows -> {draws_path}")
    print(f"Saved {len(challengers_df)} challengers -> {challengers_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
