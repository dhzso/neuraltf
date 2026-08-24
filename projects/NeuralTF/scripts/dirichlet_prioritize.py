"""Dirichlet-robust dual-track prioritization for NeuralTF candidates.

Ranks candidates by their **median integrated score across 1000 Dirichlet
weight draws** instead of the single fixed-weight score, then selects the
top 5 per track using the same composite-score bonus scheme.

SCIENTIFIC RATIONALE:
The pipeline's integrated score is a weighted sum of 7 evidence streams:
  integrated_score = sum(w_i * stream_i)
with fixed weights W = [0.211, 0.105, 0.158, 0.158, 0.105, 0.158, 0.105].

The weight choice is somewhat arbitrary. This script assesses robustness
by sampling weight vectors from a Dirichlet distribution centered on W
(concentration k=40, equivalent to ~40 pseudo-observations, giving
~95% of weight mass within ±0.1 of default weights). For each of 1000
draws, all 99 candidates are re-scored with the sampled weights. The
median score across draws is the "Dirichlet-robust" score.

This is a proper weight sensitivity analysis: the SAME weight vector is
applied to ALL candidates per draw (not per-candidate weights). NaN
streams are zeroed out in the dot product, so missing evidence contributes
zero regardless of the sampled weight.

Outputs (CSVs/MD into `projects/NeuralTF/results/`, gitignored):
  - dirichlet_top10_prioritized.csv   (5 Track A + 5 Track B)
  - dirichlet_overall_top10.csv        (overall top-10 by Dirichlet median)
  - dirichlet_overall_top10_byscore.csv (overall top-10 by score, reference)
  - dirichlet_candidate_summary_report.md

Usage:
    python projects/NeuralTF/scripts/dirichlet_prioritize.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from bioforge.projects.neuraltf.prioritize import (
    attach_v4,
    assign_tracks,
    compute_composite,
    extract_gene_symbol,
    map_v6_to_v4,
    merge_annotations,
    prepare_candidates,
    rnai_marker_notes,
    select_top,
    summarize_annotations,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO   = Path(__file__).resolve().parents[3]
DATA   = REPO / "projects" / "NeuralTF" / "data"
RUN    = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
OUT    = REPO / "projects" / "NeuralTF" / "results"

KING_DIR = REPO / "datasets" / "raw" / "Supplementary_Data_ King_2024"

# ---------------------------------------------------------------------------
# Constants (same as dirichlet_rank_analysis.py)
# ---------------------------------------------------------------------------
STREAMS = ["expression", "specificity", "reproducibility", "rnai",
           "correlation", "neural_enriched", "neural_specificity"]
W_DEFAULT = np.array([0.211, 0.105, 0.158, 0.158, 0.105, 0.158, 0.105])

N_DRAWS = 1000
K_DIR   = 40.0
SEED    = 2024


# ---------------------------------------------------------------------------
# Helpers
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
        raise ValueError(f"mmc4 header not found in {path}")
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


def integrated_scores(S: np.ndarray, W: np.ndarray) -> np.ndarray:
    """
    Compute integrated scores with Dirichlet-sampled weights, renormalizing
    over available streams for each candidate.

    This preserves the fixed-weight scorer's missing-data philosophy:
    weights are renormalized over streams available for each candidate,
    so missing evidence does not consume weight mass.

    Parameters
    ----------
    S : np.ndarray of shape (n_candidates, n_streams)
        Stream scores, with NaN for missing evidence.
    W : np.ndarray of shape (n_streams,) or (n_draws, n_streams)
        Weight vector(s). If 2D, assumed to be (n_draws, n_streams).

    Returns
    -------
    np.ndarray of shape (n_candidates,) or (n_draws, n_candidates)
        Integrated scores with missing-data renormalization.
    """
    S = np.asarray(S, dtype=float)
    W = np.asarray(W, dtype=float)

    # Handle both single weight vector (1D) and multiple draws (2D)
    if W.ndim == 1:
        W = W[np.newaxis, :]  # shape (1, n_streams)

    n_draws, n_streams = W.shape
    n_candidates = S.shape[0]

    # Mask of available evidence per candidate
    mask = ~np.isnan(S)  # shape (n_candidates, n_streams)

    # Numerator: sum over available streams of (score * weight)
    # For missing streams, treat score as 0 (so they don't contribute)
    S_filled = np.where(np.isnan(S), 0.0, S)  # shape (n_candidates, n_streams)
    numerator = np.sum(S_filled[:, np.newaxis, :] * W, axis=2)  # shape (n_candidates, n_draws)

    # Denominator: sum of weights for available streams per candidate
    # W shape: (n_draws, n_streams), mask: (n_candidates, n_streams)
    # We need to sum weights for available streams per candidate per draw
    denominator = np.sum(mask[:, np.newaxis, :] * W, axis=2)  # shape (n_candidates, n_draws)

    # Safe division: where denominator is 0, return 0
    with np.errstate(divide='ignore', invalid='ignore'):
        scores = np.where(denominator > 0, numerator / denominator, 0.0)

    return scores.T  # shape (n_draws, n_candidates) -> transpose for (n_candidates, n_draws)


def dirichlet_median_scores(S: np.ndarray, W: np.ndarray,
                            n_draws: int, k: float,
                            rng: np.random.Generator) -> np.ndarray:
    """Return the median integrated score for each candidate across draws.

    Samples ONE weight vector per draw from Dirichlet(alpha = W * k),
    applies it to ALL candidates. Missing evidence is handled by
    renormalizing weights over available streams per candidate,
    preserving the fixed-weight scorer's missing-data philosophy.
    """
    n_candidates = S.shape[0]
    all_scores = np.empty((n_draws, n_candidates), dtype=np.float32)
    mask = ~np.isnan(S)

    for d in range(n_draws):
        # Sample ONE weight vector for this draw
        alpha = W * k + 1e-9
        w = rng.gamma(W * k + 1e-9, 1.0)
        w = w / w.sum()

        # Compute integrated scores with renormalization for this draw
        mask = ~np.isnan(S)
        S_filled = np.where(np.isnan(S), 0.0, S)
        numerator = np.sum(S_filled * w, axis=1)
        denominator = np.sum(mask * w, axis=1)
        with np.errstate(divide='ignore', invalid='ignore'):
            scores = np.where(denominator > 0, numerator / denominator, 0.0)
        all_scores[d] = scores

    return np.median(all_scores, axis=0)


def compute_x1_dynamics(plass_path: Path, genes: list[str],
                        out_cache: Path) -> dict[str, float]:
    """X1 neoblast mean expression per gene, clustering on the fly."""
    import hashlib
    import scanpy as sc
    from scipy import sparse

    fingerprint = hashlib.md5(
        "\n".join(sorted(genes)).encode("utf-8")).hexdigest()[:10]
    if out_cache.exists():
        c = pd.read_csv(out_cache)
        if "fingerprint" in c.columns and str(c["fingerprint"].iloc[0]) == fingerprint:
            return dict(zip(c["gene_id"], c["mean_x1"]))

    adata = sc.read_h5ad(plass_path)
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)

    slot_of = {}
    for v in adata.var_names:
        s = str(v)
        if s.startswith("dd_Smed_v6_"):
            parts = s.split("_")
            if len(parts) >= 4:
                slot_of.setdefault(parts[3], s)
    gene_loc = {}
    for g in genes:
        num = g.split("_")[3]
        gene_loc[g] = slot_of.get(num)

    wi = slot_of.get("11973")  # smedwi-1
    sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat",
                                inplace=True)
    sc.pp.pca(adata, n_comps=50, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(adata, resolution=0.5, flavor="igraph")

    groups = pd.Series(adata.obs["leiden"].astype(str), index=adata.obs.index)
    smedwi = np.zeros(adata.n_obs)
    if wi is not None:
        smedwi = np.asarray(adata[:, wi].X.todense()).ravel()
    per = pd.Series(smedwi, index=groups.index).groupby(groups).mean()
    x1_group = per.idxmax()
    mask = (groups.values == x1_group)

    out: dict[str, float] = {}
    for g in genes:
        loc = gene_loc.get(g)
        if loc is None:
            out[g] = np.nan
            continue
        vals = np.asarray(adata[mask, loc].X.todense()).ravel()
        out[g] = float(vals.mean())
    pd.DataFrame({"gene_id": list(out), "mean_x1": list(out.values()),
                  "fingerprint": fingerprint}).to_csv(out_cache, index=False)
    return out


# ---------------------------------------------------------------------------
# Report builder (adapted from prioritize_neural_tfs.py)
# ---------------------------------------------------------------------------
def clean_ortholog(value: str, planmine_desc: str = "") -> str:
    v = str(value or "").strip()
    if v.lower() in ("nan", "none", ""):
        v = ""
    if not v:
        v = extract_gene_symbol(planmine_desc)
    return v if v.lower() not in ("nan", "none") else ""


def build_csv(top: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame()
    out["gene_id_v6"] = top["gene_id"]
    out["gene_id_v4"] = top["gene_id_v4"]
    out["gene_name"] = top["gene_name"]
    out["track"] = top["track"]
    out["rank"] = top["rank"]
    out["dirichlet_median_score"] = top["dirichlet_median_score"]
    out["composite_score"] = top["composite_score"]
    out["proof_status"] = top["proof_status"]
    out["interpro_domains"] = top["domains_all"]
    out["human_ortholog"] = [
        clean_ortholog(r["human_ortholog"], r.get("planmine_human_ortholog_desc", ""))
        for _, r in top.iterrows()
    ]
    out["rnai_phenotype_notes"] = top["rnai_phenotype_notes"]
    return out


def build_report(top: pd.DataFrame, baseline_top: pd.DataFrame,
                 g0: dict, x1: dict, n_candidates: int) -> str:
    lines: list[str] = []
    lines.append("# Dirichlet-Robust Candidate Summary Report — NeuralTF\n")
    lines.append(
        f"Inputs: `{n_candidates}` neural candidates ranked by **median "
        "integrated score** across 1000 Dirichlet weight draws "
        f"(k={K_DIR:.0f}, seed={SEED}), PlanMine annotations, King 2024 "
        "supplementary tables, v6→v4 bridge, G0 atlas.\n")
    lines.append("## Method\n")
    lines.append(
        "Same two-track scheme as the fixed-weight pipeline:\n\n"
        "- **Track A** — RNAi-validated benchmark TFs. Top 5 by composite.\n"
        "- **Track B** — novel candidates with a TF domain hit. Top 5.\n\n"
        "The key difference: the base score is the **Dirichlet median** "
        "(robust to plausible weight perturbations) rather than a single "
        "fixed-weight composite.\n")
    lines.append("## Shortlist\n")
    lines.append("| v6 id | gene_name | track | Dirichlet median | composite | human ortholog |")
    lines.append("|---|---|---|---|---|---|")
    for _, r in top.iterrows():
        ho = clean_ortholog(r["human_ortholog"], r.get("planmine_human_ortholog_desc", ""))
        lines.append(
            f"| {r['gene_id']} | {r['gene_name']} | {r['track']} "
            f"| {r['dirichlet_median_score']:.4f} "
            f"| {r['composite_score']:.3f} | {ho} |")
    lines.append("")

    # Comparison with fixed-weight baseline
    lines.append("## Comparison with fixed-weight baseline\n")
    lines.append("| gene_name | track | fixed-weight rank | Dirichlet rank | shift |")
    lines.append("|---|---|---|---|---|")
    baseline_rank_map = {}
    if baseline_top is not None and not baseline_top.empty:
        for _, r in baseline_top.iterrows():
            baseline_rank_map[r["gene_id"]] = int(r["rank"])
    for _, r in top.iterrows():
        bw = baseline_rank_map.get(r["gene_id"], "—")
        dr = int(r["rank"])
        if isinstance(bw, int):
            shift = bw - dr
            arrow = "↑" if shift > 0 else ("↓" if shift < 0 else "=")
            lines.append(f"| {r['gene_name']} | {r['track']} | {bw} | {dr} | {shift:+d} {arrow} |")
        else:
            lines.append(f"| {r['gene_name']} | {r['track']} | — | {dr} | new |")
    lines.append("")

    # Per-candidate details
    for _, r in top.iterrows():
        gid = r["gene_id"]
        lines.append(f"## {r['gene_name'] or gid} (Track {r['track']}, "
                     f"rank {r['rank']})\n")
        lines.append(f"- v6 `{gid}` · v4 `{r['gene_id_v4'] or '-'}` · "
                     f"`{r['proof_status']}`")
        lines.append(f"- Dirichlet median score: `{r['dirichlet_median_score']:.4f}`")
        lines.append(f"- composite: `{r['composite_score']:.3f}` "
                     f"(pipeline integrated `{r['integrated_score']:.3f}`, "
                     f"{int(r.get('n_streams', 0))} evidence streams)")
        lines.append(f"- DNA-binding domains (PlanMine): "
                     f"`{r.get('dna_binding_domains') or 'none annotated'}`")
        go = str(r.get("go_terms", "")).strip()
        if go:
            lines.append(f"- GO terms: {go[:400]}")
        ho = clean_ortholog(r["human_ortholog"], r.get("planmine_human_ortholog_desc", ""))
        if ho:
            lines.append(f"- Human ortholog: {ho}")
        lines.append(f"- RNAi note: {r['rnai_phenotype_notes']}")
        g = g0.get(gid)
        x = x1.get(gid)
        dyn = "G0 progenitor max log2FC "
        dyn += "n/a" if g is None or (g != g) else f"`{g:.2f}`"
        if x is not None and x == x:
            dyn += f" · X1 neoblast mean (log1p CPM) `{x:.2f}`"
        else:
            dyn += " · X1 n/a"
        lines.append(f"- Cross-stage dynamics: {dyn}")
        lines.append(
            f"- Wet-lab suggestion: design dsRNA against nt 300–800 of "
            f"`{gid}` in `datasets/processed/planmine_transcripts.fasta`; "
            f"FISH probe ≈ 800 nt antisense over the CDS region.")
        lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    # --- Load inputs --------------------------------------------------------
    rank_path = RUN / "rank_neural.csv"
    bridge_path = DATA / "bridge.csv"
    ann_path = REPO / "datasets" / "processed" / "planmine_annotations.parquet"
    mmc4_path = _resolve(KING_DIR, "mmc4")
    mmc5_path = _resolve(KING_DIR, "mmc5")

    print("== Dirichlet-robust prioritization ==")
    print(f"  rank       : {rank_path}")
    print(f"  annotations: {ann_path}")

    rank = pd.read_csv(rank_path)
    bridge = pd.read_csv(bridge_path, dtype=str)
    mmc4 = read_mmc4(mmc4_path)
    mmc5 = read_mmc5(mmc5_path)
    ann = pd.read_parquet(ann_path)
    print(f"  candidates : {len(rank)}")

    # --- Build candidate frame ----------------------------------------------
    mapping = map_v6_to_v4(bridge)
    ann_sum = summarize_annotations(ann)
    cand = prepare_candidates(rank, mmc4=mmc4)
    cand = attach_v4(cand, mapping)
    cand = merge_annotations(cand, ann_sum)
    print(f"  v6->v4 mapping: {cand['v4_mapping_flag'].value_counts().to_dict()}")

    # --- Extract stream matrix and compute Dirichlet median scores ----------
    S = cand[STREAMS].to_numpy(dtype=float)
    rng = np.random.default_rng(SEED)
    median_scores = dirichlet_median_scores(S, W_DEFAULT, N_DRAWS, K_DIR, rng)
    cand["dirichlet_median_score"] = median_scores

    # Use Dirichlet median as the ranking score (replaces integrated_score)
    cand["integrated_score"] = median_scores
    cand = compute_composite(cand)

    print(f"  Dirichlet median scores: "
          f"min={median_scores.min():.4f}, max={median_scores.max():.4f}, "
          f"mean={median_scores.mean():.4f}")

    # --- Load fixed-weight baseline for comparison --------------------------
    baseline_csv = REPO / "projects" / "NeuralTF" / "results" / "top10_neural_tfs_prioritized.csv"
    baseline_top = pd.read_csv(baseline_csv) if baseline_csv.exists() else None
    if baseline_top is not None:
        baseline_top["gene_id"] = baseline_top["gene_id_v6"]
        print(f"  baseline: {len(baseline_top)} entries from {baseline_csv.name}")

    # --- Track assignment + top 5 -------------------------------------------
    a, b = assign_tracks(cand)

    # Track B: require TF domain evidence
    b_filtered = b[
        (b["dna_binding_domains"].astype(str).str.strip() != "")
        | (b["mmc4_tf_flag"].astype(str).str.upper() == "TF")
    ]
    print(f"  Track B after TF-domain filter: {len(b_filtered)}/{len(b)}")

    ta = select_top(a, 5).assign(track="A")
    tb = select_top(b_filtered, 5).assign(track="B")
    top = pd.concat([ta, tb], ignore_index=True)

    # RNAi phenotype notes
    notes = []
    for _, r in top.iterrows():
        if r["proof_status"] == "known_rnai_validated":
            notes.append(rnai_marker_notes(mmc5, r["gene_id"]))
        else:
            notes.append("Not RNAi-tested in King 2024 mmc5; novel neural-fate candidate")
    top["rnai_phenotype_notes"] = notes

    # --- Print shortlist ----------------------------------------------------
    print("\n  === Track A (RNAi-validated) top 5 ===")
    print(top[top["track"] == "A"][
        ["gene_id", "gene_name", "dirichlet_median_score", "composite_score",
         "dna_binding_domains"]].to_string(index=False))
    print("  === Track B (novel) top 5 ===")
    print(top[top["track"] == "B"][
        ["gene_id", "gene_name", "dirichlet_median_score", "composite_score",
         "dna_binding_domains"]].to_string(index=False))

    # --- Overall top-10 by Dirichlet median score ---------------------------
    # Recompute fixed-weight scores for comparison
    fixed_S = cand[STREAMS].to_numpy(dtype=float)
    fixed_mask = ~np.isnan(fixed_S)
    cand["fixed_score"] = np.where(fixed_mask, fixed_S, 0.0) @ W_DEFAULT

    print("\n  === Overall top-10 by Dirichlet median score (all 99) ===")
    overall = cand.nlargest(10, "dirichlet_median_score")
    for _, r in overall.iterrows():
        print(f"    {r['gene_name']:>8}  dirichlet_median={r['dirichlet_median_score']:.4f}  "
              f"proof={r['proof_status']}  track={r.get('track', '?')}")

    print("\n  === Fixed-weight overall top-10 (for comparison) ===")
    fixed_overall = cand.nlargest(10, "fixed_score")
    for _, r in fixed_overall.iterrows():
        print(f"    {r['gene_name']:>8}  fixed_score={r['fixed_score']:.4f}  "
              f"proof={r['proof_status']}  track={r.get('track', '?')}")

    # Compare
    dir_set = set(overall["gene_id"].tolist())
    fix_set = set(fixed_overall["gene_id"].tolist())
    print(f"\n  Overall top-10 overlap: {len(dir_set & fix_set)}/10")
    if dir_set != fix_set:
        print(f"  Fixed-only: {[cand[cand['gene_id']==g]['gene_name'].iloc[0] for g in fix_set - dir_set]}")
        print(f"  Dirichlet-only: {[cand[cand['gene_id']==g]['gene_name'].iloc[0] for g in dir_set - fix_set]}")
    else:
        print("  Top-10 identical under both methods")

    # --- Save CSVs ----------------------------------------------------------
    OUT.mkdir(parents=True, exist_ok=True)
    csv_out = build_csv(top)
    csv_path = OUT / "dirichlet_top10_prioritized.csv"
    csv_out.to_csv(csv_path, index=False)
    print(f"\n  wrote {csv_path} ({len(csv_out)} rows)")

    # Track-based top-10 CSV (5 Track A + 5 Track B) — this is the "overall" shortlist
    top_csv = pd.DataFrame({
        "gene_id_v6": top["gene_id"],
        "gene_id_v4": top["gene_id_v4"],
        "gene_name": top["gene_name"],
        "track": top["track"],
        "rank": top["rank"],
        "dirichlet_median_score": top["dirichlet_median_score"],
        "composite_score": top["composite_score"],
        "proof_status": top["proof_status"],
    })
    top_csv_path = OUT / "dirichlet_overall_top10.csv"
    top_csv.to_csv(top_csv_path, index=False)
    print(f"  wrote {top_csv_path} (track-based 5A+5B)")

    # Actual overall top-10 by Dirichlet median score (all 99) — for reference
    overall_ref_csv = pd.DataFrame({
        "gene_id_v6": overall["gene_id"],
        "gene_name": overall["gene_name"],
        "dirichlet_median_score": overall["dirichlet_median_score"],
        "fixed_weight_score": overall["fixed_score"],
        "proof_status": overall["proof_status"],
        "track": overall.get("track", ""),
    })
    overall_ref_path = OUT / "dirichlet_overall_top10_byscore.csv"
    overall_ref_csv.to_csv(overall_ref_path, index=False)
    print(f"  wrote {overall_ref_path} (overall top-10 by score)")

    # --- Full-rank CSV for all 99 candidates ----------------------------------
    full_rank_csv = pd.DataFrame({
        "gene_id_v6": cand["gene_id"],
        "gene_name": cand["gene_name"],
        "dirichlet_median_score": cand["dirichlet_median_score"],
        "fixed_weight_score": cand["fixed_score"],
        "proof_status": cand["proof_status"],
    })
    full_rank_path = OUT / "dirichlet_centered_full_rank.csv"
    full_rank_csv.to_csv(full_rank_path, index=False)
    print(f"  wrote {full_rank_path} ({len(full_rank_csv)} rows — all 99 candidates)")

    # --- X1 dynamics --------------------------------------------------------
    x1: dict[str, float] = {}
    plass_path = REPO / "datasets" / "processed" / "plass_v6.h5ad"
    if plass_path.exists():
        try:
            x1 = compute_x1_dynamics(
                plass_path, top["gene_id"].tolist(),
                REPO / "datasets" / "processed" / "plass_x1_summary.csv")
        except Exception as exc:
            print(f"  [warn] X1 clustering failed: {exc}")

    # --- G0 atlas -----------------------------------------------------------
    g0: dict[str, float] = {}
    king_atlas = DATA / "king_atlas.tsv"
    if king_atlas.exists():
        atlas = pd.read_csv(king_atlas, sep="\t")
        g0 = atlas.groupby("v6_id")["log2fc"].max().to_dict()

    # --- Build report -------------------------------------------------------
    report = build_report(top, baseline_top, g0, x1, n_candidates=len(rank))
    report_path = OUT / "dirichlet_candidate_summary_report.md"
    report_path.write_text(report, encoding="utf-8")
    print(f"  wrote {report_path}")
    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
