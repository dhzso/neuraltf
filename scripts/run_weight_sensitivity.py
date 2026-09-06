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
  the vacated SLOT within the SAME TRACK (per-slot accounting: a Track-B
  entrant can never be credited with displacing a Track-A gene). The
  2026-09-06 audit found the previous "track-aware" code read a `track`
  column that assign_tracks() never emitted — short.get("track","?")
  returned the scalar "?" so every dropout was classified Track A and the
  displacer was the last row of the concatenated A+B shortlist; the
  shipped CSV contained impossible Track-B-displaces-Track-A attributions.
- Rank distributions: every gene that ever enters a shortlist or holds a
  top-30 rank has its FULL 1000-draw rank vector persisted (previously
  rows were kept only for in-shortlist/top-30 draws — a truncation that
  biased fig 06's stability boxplots, e.g. Zeb-1 median 49 vs true 76).
- Determinism: boundary ties are resolved by the same tie-break column
  set used by select_top (draw base score -> integrated -> n_streams ->
  gene_id ascending), never by unstable argsort artifact.

Draw weights come from a symmetric Dirichlet(1,...,1) (uniform over the
weight simplex); per-draw scores renormalize over each candidate's
available streams, mirroring EvidenceScorer exactly.

Outputs (written to projects/NeuralTF/figures/ for figures 06/07):
- weight_sensitivity_draws.csv: draw, gene_id, gene_name, rank, in_top_10
  (complete 1000-draw rank vectors for entrants/top-30 genes)
- weight_sensitivity_top10_challengers.csv: per-candidate summary
  (baseline_rank, frac_draws_in_top10, best/median rank over ALL draws,
  per-slot same-track displacement info)
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
# Tie-breaks mirror select_top: the method's OWN base score first (for a
# draw, _draw_base), then integrated_score, then n_streams, then gene_id
# ascending. The previous TIE_COLS omitted the draw base, so exact
# composite ties could resolve differently here than in the published
# selection.
TIE_COLS = ["_draw_base", "integrated_score", "n_streams"]


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
    Track A/B split, Track-B domain gate, deterministic tie-breaks.

    The returned frame carries a real `track` column ("A"/"B") and a
    `rank_within_track` (1..5) — required by the per-slot displacement
    attribution. assign_tracks() filters by proof_status but never
    annotates the track, which is why this function does it here.
    """
    d = cand.copy()
    d["_draw_base"] = scores
    d = apply_bonuses(d, "_draw_base")
    a, b = assign_tracks(d)
    b = gate_track_b(b)

    def top5(track: pd.DataFrame, label: str) -> pd.DataFrame:
        if track.empty:
            return track
        tmp = track.copy()
        tmp["_gid_desc"] = [_gene_id_desc(g) for g in tmp["gene_id"]]
        cols = ["composite_score"] + [c for c in TIE_COLS if c in tmp.columns] + ["_gid_desc"]
        out = tmp.sort_values(cols, ascending=False).head(5)
        out["track"] = label
        out["rank_within_track"] = range(1, len(out) + 1)
        return out

    return pd.concat([top5(a, "A"), top5(b, "B")], ignore_index=True)


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

    # Baseline shortlist with REAL track labels (from the published CSV,
    # which carries a `track` column; recomputed shortlists get theirs
    # from build_shortlist).
    if baseline_path.exists():
        baseline_df = pd.read_csv(baseline_path)
        gid_col = "gene_id_v6" if "gene_id_v6" in baseline_df.columns else "gene_id"
        baseline_shortlist = set(baseline_df[gid_col].astype(str))
        baseline_track = dict(zip(
            baseline_df[gid_col].astype(str),
            baseline_df["track"].astype(str) if "track" in baseline_df.columns
            else ["?"] * len(baseline_df)))
    else:
        base_short = build_shortlist(
            cand, weighted_scores(S, np.full(len(available), 1 / len(available))))
        baseline_shortlist = set(base_short["gene_id"])
        baseline_track = dict(zip(base_short["gene_id"], base_short["track"]))
    print(f"Baseline shortlist: {len(baseline_shortlist)} genes")

    baseline_rank_map = dict(zip(
        cand["gene_id"],
        pd.to_numeric(cand["integrated_score"], errors="coerce").fillna(0)
        .rank(ascending=False, method="min").astype(int)))

    draw_rows = []
    in_top10_count = {g: 0 for g in cand["gene_id"]}
    best_rank = {g: n + 1 for g in cand["gene_id"]}
    # Full per-draw rank vectors for EVERY gene (held in memory; the
    # persisted CSV subset is materialized AFTER the loop from these
    # complete vectors, so a gene that first hits the top-30 at draw 700
    # still gets all 1000 rows — a mid-run persistence flag would
    # silently truncate its early draws).
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
        # Per-slot, same-track challenger attribution. A baseline gene
        # that misses the draw's shortlist is credited against the
        # draw-NEW gene that entered ITS OWN track (matched by slot
        # order when the entrant count allows; otherwise the entrant's
        # lowest vacated slot). Cross-track attribution is structurally
        # impossible: entrants and dropouts are paired within track.
        draw_track_of = dict(zip(short["gene_id"], short["track"]))
        entrants_by_track: dict[str, list[str]] = {"A": [], "B": []}
        for _, srow in short.iterrows():
            gid = str(srow["gene_id"])
            if gid not in baseline_shortlist:
                entrants_by_track[str(srow["track"])].append(gid)
        dropped_by_track: dict[str, list[str]] = {"A": [], "B": []}
        for gid in baseline_shortlist - short_set:
            bt = baseline_track.get(gid, "?")
            dropped_by_track[bt if bt in ("A", "B") else "A"].append(gid)

        # Slot-order pairing within each track: entrants sorted by their
        # rank_within_track in this draw fill the vacated slots of the
        # dropped genes ordered by their baseline rank within track.
        base_rank_in_track: dict[str, dict[str, int]] = {"A": {}, "B": {}}
        if baseline_path.exists():
            for _, brow in baseline_df.iterrows():
                g = str(brow[gid_col])
                t = str(brow.get("track", "?"))
                rk = brow.get("rank", 99)
                if t in base_rank_in_track and g in baseline_shortlist:
                    try:
                        base_rank_in_track[t][g] = int(rk)
                    except (TypeError, ValueError):
                        base_rank_in_track[t][g] = 99
        else:
            for gid, t in baseline_track.items():
                base_rank_in_track[t][gid] = 99

        attribution: dict[str, str] = {}
        for t in ("A", "B"):
            drops = sorted(dropped_by_track[t],
                           key=lambda g: base_rank_in_track[t].get(g, 99))
            ents = [g for g in entrants_by_track[t]]
            ent_rank = {g: int(short.loc[short["gene_id"] == g,
                                          "rank_within_track"].iloc[0])
                        for g in ents}
            ents = sorted(ents, key=lambda g: ent_rank.get(g, 99))
            for drop_g, ent_g in zip(drops, ents):
                attribution[drop_g] = ent_g
            # More dropouts than entrants: remaining dropouts are
            # attributed to the last entrant of the same track (the gene
            # holding the marginal seat) — never to the other track.
            for extra_drop in drops[len(ents):]:
                if ents:
                    attribution[extra_drop] = ents[-1]

        for i, gid in enumerate(cand["gene_id"]):
            in10 = gid in short_set
            rank_samples[gid].append(int(ranks[i]))
            if int(ranks[i]) < best_rank[gid]:
                best_rank[gid] = int(ranks[i])
            if in10:
                in_top10_count[gid] += 1
            elif gid in attribution:
                displacer = attribution[gid]
                displaced_by[gid][displacer] = displaced_by[gid].get(displacer, 0) + 1

        if (d + 1) % 100 == 0:
            print(f"  {d + 1}/{N_DRAWS}", flush=True)

    # Materialize the draws CSV from the COMPLETE in-memory rank vectors:
    # every gene that ever entered a shortlist or held a top-30 rank gets
    # its full 1000-row vector (persistence decided after the loop, when
    # all eligibility information is available — no mid-run truncation).
    keep_genes = {g for g in cand["gene_id"]
                  if in_top10_count[g] > 0 or (rank_samples[g] and min(rank_samples[g]) <= 30)}
    gene_names = {g: cand["gene_name"].iloc[i]
                  for i, g in enumerate(cand["gene_id"])}
    top30_overall = {g for g in cand["gene_id"]
                     if rank_samples[g] and min(rank_samples[g]) <= 30}
    shortlist_ever = {g for g in cand["gene_id"] if in_top10_count[g] > 0}
    draw_rows = []
    for gid in sorted(keep_genes):
        name = gene_names.get(gid, "")
        for d_i, rk in enumerate(rank_samples[gid]):
            draw_rows.append((d_i + 1, gid, name, int(rk),
                              gid in shortlist_ever))

    # Cross-track attribution guard: after the rewrite it must be
    # structurally impossible for a displacer to sit in a different
    # track than the gene it displaced.
    for gid, disp in [(g, d) for g, dd in displaced_by.items() for d in dd]:
        gt = baseline_track.get(gid, "?")
        dt = baseline_track.get(disp, "?")
        if gt in ("A", "B") and dt in ("A", "B") and gt != dt:
            raise AssertionError(
                f"cross-track attribution leaked through: {gid} ({gt}) "
                f"attributed to {disp} ({dt})")

        if (d + 1) % 100 == 0:
            print(f"  {d + 1}/{N_DRAWS}", flush=True)

    draws_df = pd.DataFrame(draw_rows, columns=["draw", "gene_id", "gene_name", "rank", "in_top_10"])

    # Completeness guard: every persisted gene must have exactly N_DRAWS
    # rows (full rank vector). A partial vector would silently re-introduce
    # the truncation bias in any downstream boxplot.
    _counts = draws_df.groupby("gene_id").size()
    _incomplete = _counts[_counts != N_DRAWS]
    if len(_incomplete) > 0:
        raise AssertionError(
            f"truncated rank vectors detected for {len(_incomplete)} genes "
            f"(e.g. {_incomplete.head(3).to_dict()}) - persistence rule is broken")

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
            "baseline_track": baseline_track.get(gid, ""),
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
