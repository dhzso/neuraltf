#!/usr/bin/env python
"""Export ranked TF (Transcription Factor) CSVs.

Three scope levels from King 2024 mmc4 TF catalog:

  19 TFs — neural-filtered: TFs with King neural signal or RNAi evidence
  43 TFs — candidates: TFs that passed expression filter (p ≤ 0.05)
  74 TFs — catalog: all TFs from King mmc4 TF sheet

All outputs sorted by composite score (descending), with rich columns:
  gene_id_v6, gene_id_v4, gene_name, track, rank, composite_score,
  proof_status, interpro_domains, human_ortholog, rnai_phenotype_notes

Outputs:
  projects/NeuralTF/results/tf_ranked_neural_top19.csv
  projects/NeuralTF/results/tf_ranked_all_top43.csv
  projects/NeuralTF/results/tf_ranked_catalog_top74.csv

Usage:
    python projects/NeuralTF/scripts/export_fstf_ranked.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / "projects" / "NeuralTF" / "data"
RUN  = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
OUT  = REPO / "projects" / "NeuralTF" / "results"
OUT.mkdir(parents=True, exist_ok=True)

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


def _short_dd(gene_id: str) -> str | None:
    """Extract the dd<numeric> short ID from a structured v6 ID.

    Mirrors pipeline._short_id: 'dd_Smed_v6_11150_0_1' -> 'dd11150'.
    (mmc5's first column holds bare dd#### tokens, NOT full v6 IDs —
    the 2026-09-06 audit found the previous direct-equality check
    `gene_id not in mmc5['fstf_rnai'].values` was ALWAYS true, so every
    gene, including RNAi-validated ones, exported as "Not RNAi-tested".)
    """
    import re
    if not gene_id or not isinstance(gene_id, str):
        return None
    m = re.search(r"dd_Smed_v[46]_(\d+)", gene_id)
    if m:
        return f"dd{m.group(1)}"
    m = re.search(r"\bdd(\d+)\b", gene_id)
    if m:
        return f"dd{m.group(1)}"
    return None


def rnai_marker_notes(mmc5: pd.DataFrame, gene_id: str) -> str:
    if mmc5.empty:
        return "Not RNAi-tested in King 2024 mmc5"
    short = _short_dd(gene_id)
    # mmc5 col 0 stores 'dd####' (optionally with a parenthetical symbol);
    # match on the extracted short ID, never on the full v6 dialect.
    targets = mmc5["fstf_rnai"].astype(str).str.strip()
    if short:
        hit_mask = targets.apply(
            lambda t: bool(t == short
                           or (short and t.split(" (")[0].strip() == short)))
    else:
        hit_mask = pd.Series(False, index=mmc5.index)
    if not hit_mask.any():
        return "Not RNAi-tested in King 2024 mmc5"
    row = mmc5[hit_mask].iloc[0]
    markers = []
    for c in mmc5.columns[1:]:
        val = row[c]
        if pd.notna(val) and str(val).strip() not in ("", "0"):
            markers.append(f"{c.replace('marker_', '')}={val}")
    if markers:
        return f"Brain RNAi; markers: {', '.join(markers)}"
    return "Brain RNAi; no markers"


def main() -> int:
    print("== Export ranked FSTF CSVs ==")

    # --- Load data ---------------------------------------------------------
    rank = pd.read_csv(RUN / "rank.csv")
    bridge = pd.read_csv(DATA / "bridge.csv", dtype=str)
    mmc4_path = _resolve(KING_DIR, "mmc4")
    mmc4 = read_mmc4(mmc4_path) if mmc4_path.exists() else pd.DataFrame()
    mmc5_path = _resolve(KING_DIR, "mmc5")
    mmc5 = read_mmc5(mmc5_path) if mmc5_path.exists() else pd.DataFrame()
    ann_path = REPO / "datasets" / "processed" / "planmine_annotations.parquet"
    ann = pd.read_parquet(ann_path) if ann_path.exists() else pd.DataFrame()

    # --- Get FSTF IDs from mmc4 TF sheet ----------------------------------
    fstf_ids_74 = set()
    if mmc4_path.exists():
        tf_catalog = pd.read_excel(mmc4_path, sheet_name="TF")
        fstf_ids_74 = set(
            tf_catalog.loc[
                tf_catalog["FSTF?"].astype(str).str.strip().str.lower() == "yes",
                "Gene ID",
            ]
        )
    print(f"  FSTFs in catalog: {len(fstf_ids_74)}")
    print(f"  candidates: {len(rank)}")

    # --- Build candidate frame ---------------------------------------------
    from bioforge.projects.neuraltf.prioritize import (
        map_v6_to_v4, prepare_candidates, attach_v4, merge_annotations,
        summarize_annotations, assign_tracks, compute_composite,
    )

    mapping = map_v6_to_v4(bridge)
    ann_sum = summarize_annotations(ann) if not ann.empty else pd.DataFrame()
    cand = prepare_candidates(rank, mmc4=mmc4)
    cand = attach_v4(cand, mapping)
    if not ann_sum.empty:
        cand = merge_annotations(cand, ann_sum)
    cand = compute_composite(cand)

    # --- Track assignment ---------------------------------------------------
    a, b = assign_tracks(cand)
    b_tf = b[
        (b["dna_binding_domains"].astype(str).str.strip() != "")
        | (b["mmc4_tf_flag"].astype(str).str.upper() == "TF")
    ]
    cand.loc[a.index, "track"] = "A"
    cand.loc[b_tf.index, "track"] = "B"
    other = cand[cand["track"].isna()]
    cand.loc[other.index, "track"] = "-"

    cand = cand.sort_values("composite_score", ascending=False).reset_index(drop=True)
    cand["rank"] = range(1, len(cand) + 1)

    # --- RNAi phenotype notes -----------------------------------------------
    notes = []
    for _, r in cand.iterrows():
        if r["proof_status"] == "known_rnai_validated":
            notes.append(rnai_marker_notes(mmc5, r["gene_id"]))
        else:
            notes.append("Not RNAi-tested in King 2024 mmc5; novel candidate")
    cand["rnai_phenotype_notes"] = notes

    # Rename the underlying column to match the harmonized header name
    # (the pipeline stores domains under 'domains_all'; the publication
    # header is 'interpro_domains' for descriptive accuracy).
    if "domains_all" in cand.columns and "interpro_domains" not in cand.columns:
        cand = cand.rename(columns={"domains_all": "interpro_domains"})

    out_cols = [
        "gene_id_v6", "gene_id_v4", "gene_name", "track", "rank",
        "composite_score", "proof_status", "interpro_domains",
        "human_ortholog", "rnai_phenotype_notes",
    ]

    # --- Scope 1: 19 TFs in neural-filtered set --------------------------
    neural_mask = cand["neural_enriched"].notna() | (cand["rnai"] > 0)
    fstf_19 = cand[neural_mask & cand["gene_id"].isin(fstf_ids_74)].copy()
    fstf_19["rank"] = range(1, len(fstf_19) + 1)
    fstf_19[out_cols].to_csv(OUT / "tf_ranked_neural_top19.csv", index=False)
    print(f"  wrote tf_ranked_neural_top19.csv ({len(fstf_19)} rows)")

    # --- Scope 2: 43 TFs in all candidates -------------------------------
    fstf_43 = cand[cand["gene_id"].isin(fstf_ids_74)].copy()
    fstf_43["rank"] = range(1, len(fstf_43) + 1)
    fstf_43[out_cols].to_csv(OUT / "tf_ranked_all_top43.csv", index=False)
    print(f"  wrote tf_ranked_all_top43.csv ({len(fstf_43)} rows)")

    # --- Scope 3: 74 TFs from catalog (full list) ------------------------
    # For catalog FSTFs not in 249 candidates, create rows with available data
    fstf_74_in_cand = cand[cand["gene_id"].isin(fstf_ids_74)].copy()
    fstf_74_in_cand["rank"] = range(1, len(fstf_74_in_cand) + 1)

    # Catalog FSTFs not in candidates (no expression data)
    fstf_74_missing_ids = fstf_ids_74 - set(fstf_74_in_cand["gene_id"])
    if fstf_74_missing_ids:
        missing_rows = []
        for gid in sorted(fstf_74_missing_ids):
            missing_rows.append({
                "gene_id_v6": gid,
                "gene_id_v4": "",
                "gene_name": gid.replace("dd_Smed_v6_", "dd").replace("_0_1", "").replace("_1_1", ""),
                "track": "-",
                "rank": len(fstf_74_in_cand) + len(missing_rows) + 1,
                "composite_score": float("nan"),
                "proof_status": "catalog_fstf_not_in_candidates",
                "interpro_domains": "",
                "human_ortholog": "",
                "rnai_phenotype_notes": "",
            })
        fstf_74_missing = pd.DataFrame(missing_rows)
        fstf_74 = pd.concat([fstf_74_in_cand[out_cols], fstf_74_missing], ignore_index=True)
    else:
        fstf_74 = fstf_74_in_cand[out_cols]

    fstf_74.to_csv(OUT / "tf_ranked_catalog_top74.csv", index=False)
    print(f"  wrote tf_ranked_catalog_top74.csv ({len(fstf_74)} rows)")

    print(f"\n  Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
