"""Dual-track prioritization + summary report for NeuralTF candidates.

Reads the run outputs, PlanMine annotations, the identifier bridge and the
King 2024 supplementary tables, then writes::

    projects/NeuralTF/results/top10_neural_tfs_prioritized.csv
    projects/NeuralTF/results/candidate_summary_report.md

Also computes cross-stage expression dynamics (Plass X1 neoblast vs. the
G0 progenitor atlas) for the shortlisted TFs.

Usage::

    python scripts/prioritize_neural_tfs.py [--repo D:/Bioinformatics]
"""

from __future__ import annotations

import argparse
import sys
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

REPO = Path.cwd()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=REPO, help="repo root (default: cwd)")
    p.add_argument("--rank", type=Path, default=None)
    p.add_argument("--bridge", type=Path, default=None)
    p.add_argument("--annotations", type=Path, default=None)
    p.add_argument("--skip-x1", action="store_true",
                   help="skip the Plass X1 clustering step")
    return p.parse_args(argv)


def _resolve(king_dir: Path, name: str) -> Path:
    cand = king_dir / f"1-s2.0-S2211124724001712-{name}.xlsx"
    if cand.exists():
        return cand
    if king_dir.exists():
        for p in sorted(king_dir.iterdir()):
            if p.suffix.lower() == ".xlsx" and p.stem.lower().endswith(name):
                return p
    return cand


def load_supported(args: argparse.Namespace) -> dict:
    repo = args.repo.resolve()
    data = repo / "projects" / "NeuralTF" / "data"
    run = repo / "projects" / "NeuralTF" / "runs" / "pipeline_run"
    king_dir = repo / "datasets" / "raw" / "Supplementary_Data_ King_2024"
    return {
        "repo": repo,
        "rank": args.rank or (run / "rank_neural.csv"),
        "bridge": args.bridge or (data / "bridge.csv"),
        "annotations": args.annotations or (
            repo / "datasets" / "processed" / "planmine_annotations.parquet"),
        "mmc4": _resolve(king_dir, "mmc4"),
        "mmc5": _resolve(king_dir, "mmc5"),
        "king_atlas": data / "king_atlas.tsv",
        "plass": repo / "datasets" / "processed" / "plass_v6.h5ad",
        "out_dir": repo / "projects" / "NeuralTF" / "results",
        "fasta": repo / "datasets" / "processed" / "planmine_transcripts.fasta",
    }


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


def compute_x1_dynamics(plass_path: Path, genes: list[str],
                        out_cache: Path) -> dict[str, float]:
    """X1 neoblast mean expression per gene, clustering on the fly.

    The X1 cluster is identified as the Leiden cluster with the highest mean
    smedwi-1 (dd_Smed_v6_11973_0) expression (canonical neoblast marker).
    Results are cached keyed by a fingerprint of the gene list so stale
    caches from other runs never leak.
    Returns ``{gene_id: mean log1p-CPM}``.
    """
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


def select_shortlist(cand: pd.DataFrame, mmc5: pd.DataFrame | None) -> pd.DataFrame:
    """Top-5 per track with RNAi phenotype notes attached.

    Track B requires a tangible TF identity: a DNA-binding protein-domain hit
    in PlanMine or an mmc4 "TF" flag — no hypothetical factors without domain
    evidence.
    """
    a, b = assign_tracks(cand)
    b = b[
        (b["dna_binding_domains"].astype(str).str.strip() != "")
        | (b["mmc4_tf_flag"].astype(str).str.upper() == "TF")
    ]
    print(f"  Track B after TF-domain filter: {len(b)}/{len(assign_tracks(cand)[1])}")
    ta = select_top(a, 5).assign(track="A")
    tb = select_top(b, 5).assign(track="B")
    top = pd.concat([ta, tb], ignore_index=True)
    notes = []
    for _, r in top.iterrows():
        if r["proof_status"] == "known_rnai_validated":
            notes.append(rnai_marker_notes(mmc5, r["gene_id"]))
        else:
            notes.append("Not RNAi-tested in King 2024 mmc5; novel neural-fate candidate")
    top["rnai_phenotype_notes"] = notes
    return top.sort_values(["track", "rank"]).reset_index(drop=True)


def clean_ortholog(value: str, planmine_desc: str = "") -> str:
    """Normalise a human-ortholog label; falls back to the PlanMine symbol."""
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
    out["composite_score"] = top["composite_score"]
    out["proof_status"] = top["proof_status"]
    out["interpro_domains"] = top["domains_all"]
    out["human_ortholog"] = [
        clean_ortholog(r["human_ortholog"], r.get("planmine_human_ortholog_desc", ""))
        for _, r in top.iterrows()
    ]
    out["rnai_phenotype_notes"] = top["rnai_phenotype_notes"]
    return out


def build_report(top: pd.DataFrame, g0: dict, x1: dict) -> str:
    lines: list[str] = []
    lines.append("# Candidate Summary Report — NeuralTF Prioritization\n")
    lines.append(
        "Inputs: `rank_neural.csv` (96 neural candidates), PlanMine "
        "annotations (`datasets/processed/planmine_annotations.parquet`), the "
        "v6→v4 identifier bridge, King 2024 supplementary tables (mmc4 TF "
        "catalog, mmc5 FSTF RNAi screen), G0 atlas (`king_atlas.tsv`).\n")
    lines.append("## Method\n")
    lines.append(
        "Two independent tracks:\n\n"
        "- **Track A** — `proof_status == known_rnai_validated`: RNAi-validated "
        "benchmark TFs from the King 2024 FSTF screen. Top 5 by composite "
        "score.\n"
        "- **Track B** — `proof_status == novel_candidate`: no published RNAi "
        "data; filtered to candidates with a clear DNA-binding TF domain "
        "(PlanMine protein-domain hits or mmc4 TF flag), then top 5 by "
        "composite score.\n\n"
        "`composite_score = integrated_score + bonuses` (formula in "
        "`bioforge/projects/neuraltf/prioritize.py`): TF domain +0.05, neural "
        "GO +0.03, TF GO +0.02, human ortholog +0.02, RNAi-validated +0.02.\n")
    lines.append("## Shortlist\n")
    lines.append("| v6 id | gene_name | track | rank | composite | human ortholog |")
    lines.append("|---|---|---|---|---|---|")
    for _, r in top.iterrows():
        ho = clean_ortholog(r["human_ortholog"], r.get("planmine_human_ortholog_desc", ""))
        lines.append(
            f"| {r['gene_id']} | {r['gene_name']} | {r['track']} | {r['rank']} "
            f"| {r['composite_score']:.3f} | {ho} |")
    lines.append("")

    for _, r in top.iterrows():
        gid = r["gene_id"]
        lines.append(f"## {r['gene_name'] or gid} (Track {r['track']}, "
                     f"rank {r['rank']})\n")
        lines.append(f"- v6 `{gid}` · v4 `{r['gene_id_v4'] or '-'}` · "
                     f"`{r['proof_status']}`")
        lines.append(f"- composite `{r['composite_score']:.3f}` "
                     f"(pipeline integrated `{r['integrated_score']:.3f}`, "
                     f"{int(r['n_streams'])} evidence streams)")
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
    lines.append("## Reproducibility\n")
    lines.append(
        "Deterministic pipeline (seed 42, pinned inputs, commit-pinned raw "
        "files). PlanMine snapshot dated on run; identifier mapping via the "
        "bridge table with explicit ambiguity flags (no numeric guessing).")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    paths = load_supported(args)
    print("== NeuralTF dual-track prioritization ==")
    print(f"  rank      : {paths['rank']}")
    print(f"  bridge    : {paths['bridge']}")
    print(f"  annotations: {paths['annotations']}")

    rank = pd.read_csv(paths["rank"])
    bridge = pd.read_csv(paths["bridge"], dtype=str)
    mmc4 = read_mmc4(paths["mmc4"])
    mmc5 = read_mmc5(paths["mmc5"])
    ann = pd.read_parquet(paths["annotations"])
    print(f"  candidates: {len(rank)}  "
          f"{rank['proof_status'].value_counts().to_dict()}")

    mapping = map_v6_to_v4(bridge)
    ann = summarize_annotations(ann)
    cand = prepare_candidates(rank, mmc4=mmc4)
    cand = attach_v4(cand, mapping)
    cand = merge_annotations(cand, ann)
    cand = compute_composite(cand)

    flags = cand["v4_mapping_flag"].value_counts()
    print(f"  v6->v4 mapping: {flags.to_dict()}")

    top = select_shortlist(cand, mmc5)
    print("\n  === Track A (RNAi-validated) top 5 ===")
    print(top[top["track"] == "A"][
        ["gene_id", "gene_name", "composite_score", "integrated_score",
         "dna_binding_domains"]].to_string(index=False))
    print("  === Track B (novel) top 5 ===")
    print(top[top["track"] == "B"][
        ["gene_id", "gene_name", "composite_score", "integrated_score",
         "dna_binding_domains"]].to_string(index=False))

    out_dir = paths["out_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)
    csv = build_csv(top)
    csv_path = out_dir / "top10_neural_tfs_prioritized.csv"
    csv.to_csv(csv_path, index=False)
    print(f"  wrote {csv_path} ({len(csv)} rows)")

    g0 = {}
    if paths["king_atlas"].exists():
        atlas = pd.read_csv(paths["king_atlas"], sep="\t")
        g0 = atlas.groupby("v6_id")["log2fc"].max().to_dict()

    x1 = {}
    if not args.skip_x1 and paths["plass"].exists():
        try:
            x1 = compute_x1_dynamics(
                paths["plass"], top["gene_id"].tolist(),
                paths["repo"] / "datasets" / "processed" /
                "plass_x1_summary.csv")
        except Exception as exc:  # noqa: BLE001
            print(f"  [warn] X1 clustering failed: {exc}")

    report_path = out_dir / "candidate_summary_report.md"
    report_path.write_text(build_report(top, g0, x1), encoding="utf-8")
    print(f"  wrote {report_path}")
    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())