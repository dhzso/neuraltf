#!/usr/bin/env python
"""Weight-sensitivity analysis: 1000 random weight draws over the 99 neural TFs.

For each draw, new weights are sampled from a symmetric Dirichlet(1,...,1)
distribution (uniform over the weight simplex). Candidates are scored with
the random weights (renormalized over available streams per candidate),
ranked, and membership in the top-10 is recorded.

Outputs (written to projects/NeuralTF/figures/ so figure 06/07 can read them):
- weight_sensitivity_draws.csv: draw, gene_id, gene_name, rank, in_top_10
  (one row per draw x candidate)
- weight_sensitivity_top10_challengers.csv: per-candidate summary stats
  (baseline_rank, frac_draws_in_top10, best/median rank, displacement info)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT / "projects" / "NeuralTF" / "runs" / "pipeline_run"
FIG_DIR = ROOT / "projects" / "NeuralTF" / "figures"
RESULTS_DIR = ROOT / "projects" / "NeuralTF" / "results"

STREAMS = ["expression", "specificity", "reproducibility", "rnai",
           "correlation", "neural_enriched", "neural_specificity"]
N_DRAWS = 1000
SEED = 42


def weighted_scores(S: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Row-wise weighted scores with renormalization over available streams."""
    mask = ~np.isnan(S)
    filled = np.where(mask, S, 0.0)
    denom = (mask * w).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        num = filled @ w
        return np.where(denom > 0, num / denom, 0.0)


def main() -> int:
    neural_path = RUN_DIR / "rank_neural.csv"
    baseline_path = RESULTS_DIR / "top10_neural_tfs_prioritized.csv"
    if not neural_path.exists():
        print(f"Missing {neural_path} — run the pipeline first")
        return 1

    cand = pd.read_csv(neural_path)
    S = cand[STREAMS].to_numpy(dtype=float)
    n = len(cand)
    print(f"{n} neural candidates, {N_DRAWS} random weight draws")

    rng = np.random.default_rng(SEED)
    W = rng.dirichlet(np.ones(len(STREAMS)), size=N_DRAWS)  # (draws, 7) uniform simplex

    baseline_rank_map = dict(zip(
        cand["gene_id"], cand["integrated_score"].rank(ascending=False, method="min").astype(int)))
    baseline_top10 = set(
        pd.read_csv(baseline_path)["gene_id_v6"].tolist()) if baseline_path.exists() else set()

    draw_rows = []
    in_top10_count = {g: 0 for g in cand["gene_id"]}
    best_rank = {g: n + 1 for g in cand["gene_id"]}
    rank_samples = {g: [] for g in cand["gene_id"]}
    displaced_by = {g: {} for g in cand["gene_id"]}   # gene -> {displacer: count}

    for d in range(N_DRAWS):
        w = W[d]
        scores = weighted_scores(S, w)
        order = np.argsort(-scores)
        ranks = np.empty(n, dtype=int)
        ranks[order] = np.arange(1, n + 1)

        top10_idx = set(order[:10])
        for i, gid in enumerate(cand["gene_id"]):
            in10 = i in top10_idx
            draw_rows.append((d + 1, gid, cand["gene_name"].iloc[i], int(ranks[i]), in10))
            rank_samples[gid].append(int(ranks[i]))
            if int(ranks[i]) < best_rank[gid]:
                best_rank[gid] = int(ranks[i])
            if in10:
                in_top10_count[gid] += 1
            elif gid in baseline_top10:
                # baseline top-10 gene displaced: record which gene took the top spot this draw
                displacer = str(cand["gene_id"].iloc[order[0]])
                displaced_by[gid][displacer] = displaced_by[gid].get(displacer, 0) + 1

        if (d + 1) % 100 == 0:
            print(f"  {d + 1}/{N_DRAWS}", flush=True)

    draws_df = pd.DataFrame(draw_rows, columns=["draw", "gene_id", "gene_name", "rank", "in_top_10"])

    # Per-candidate summary: only candidates that ever entered a top-10
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
            "in_baseline_top10": gid in baseline_top10,
            "n_draws_in_top10": in_top10_count[gid],
            "frac_draws_in_top10": round(in_top10_count[gid] / N_DRAWS, 4),
            "best_rank": best_rank[gid],
            "median_rank": int(np.median(rank_samples[gid])),
            "displaced_most_often": most_displacer,
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
