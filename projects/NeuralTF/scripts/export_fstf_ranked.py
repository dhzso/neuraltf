#!/usr/bin/env python
"""Export ranked FSTF (Planarian Stem Cell Transcription Factor) CSV.

Reads the pipeline output (rank.csv - all 249 candidates) and produces
a CSV with rich columns, sorted by composite score (descending).

Outputs:
  projects/NeuralTF/results/fstf_ranked_all.csv       (all 249)
  projects/NeuralTF/results/fstf_ranked_neural.csv     (99 neural-filtered)

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


def rnai_marker_notes(mmc5: pd.DataFrame, gene_id: str) -> str:
    if mmc5.empty or gene_id not in mmc5["fstf_rnai"].values:
        return "Not RNAi-tested in King 2024 mmc5"
    row = mmc5[mmc5["fstf_rnai"] == gene_id].iloc[0]
    markers = []
    for c in mmc5.columns[1:]:
        val = row[c]
        if pd.notna(val) and str(val).strip() not in ("", "0"):
            markers.append(f"{c.replace('marker_', '')}={val}")
    if markers:
        return f"Brain RNAi; markers: {', '.join(markers)}"
    return "Brain RNAi; no markers"


def main() -> int:
    print("== Export ranked FSTF CSV ==")

    rank = pd.read_csv(RUN / "rank.csv")
    bridge = pd.read_csv(DATA / "bridge.csv", dtype=str)
    mmc4_path = _resolve(KING_DIR, "mmc4")
    mmc4 = read_mmc4(mmc4_path) if mmc4_path.exists() else pd.DataFrame()
    mmc5_path = _resolve(KING_DIR, "mmc5")
    mmc5 = read_mmc5(mmc5_path) if mmc5_path.exists() else pd.DataFrame()
    ann_path = REPO / "datasets" / "processed" / "planmine_annotations.parquet"
    ann = pd.read_parquet(ann_path) if ann_path.exists() else pd.DataFrame()

    print(f"  candidates: {len(rank)}")

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

    notes = []
    for _, r in cand.iterrows():
        if r["proof_status"] == "known_rnai_validated":
            notes.append(rnai_marker_notes(mmc5, r["gene_id"]))
        else:
            notes.append("Not RNAi-tested in King 2024 mmc5; novel candidate")
    cand["rnai_phenotype_notes"] = notes

    out_cols = [
        "gene_id_v6", "gene_id_v4", "gene_name", "track", "rank",
        "composite_score", "proof_status", "domains_all",
        "human_ortholog", "rnai_phenotype_notes",
    ]

    csv_all = cand[out_cols].copy()
    csv_all.to_csv(OUT / "fstf_ranked_all.csv", index=False)
    print(f"  wrote fstf_ranked_all.csv ({len(csv_all)} rows)")

    neural_mask = cand["neural_enriched"].notna() | (cand["rnai"] > 0)
    csv_neural = cand.loc[neural_mask, out_cols].copy()
    csv_neural["rank"] = range(1, len(csv_neural) + 1)
    csv_neural.to_csv(OUT / "fstf_ranked_neural.csv", index=False)
    print(f"  wrote fstf_ranked_neural.csv ({len(csv_neural)} rows)")

    print(f"\n  Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
