#!/usr/bin/env python
"""Build unified Master TF Catalog by merging King 2024 + Perez 2025 TF annotations.

Merges on dd_Smed_v6 gene IDs, deduplicates, flags the source of each entry
(king_only / perez_only / both), and writes a single authoritative catalog used
downstream by the NeuralTF pipeline for candidate gating and scoring.

Inputs:
  datasets/raw/Supplementary_Data_ King_2024/1-s2.0-S2211124724001712-mmc4.xlsx
  projects/NeuralTF/data/perez_tf_summary.csv   (pre-processed from MOESM5)

Output:
  projects/NeuralTF/data/master_tf_catalog.csv

Schema of output:
  v6_id          : dd_Smed_v6 gene identifier (primary key)
  gene_name      : human-readable gene name (King GenBank preferred, Perez fallback)
  tf_family_king : TF family from King mmc4 (e.g., Homeobox, bHLH)
  tf_class_perez : TF class from Perez MOESM5 (e.g., Homeodomain, bHLH)
  tf_family_perez: TF family from Perez MOESM5
  pfam_domain    : PFAM domain name (Perez)
  human_ortholog : human ortholog symbol (Perez preferred, King fallback)
  is_fstf        : True if flagged as FSTF in King mmc4
  has_rnai       : True if appears in King mmc5 RNAi screen
  in_king        : True if present in King mmc4
  in_perez       : True if present in Perez MOESM5 (with v6 mapping)
  source         : 'both' / 'king_only' / 'perez_only'

Usage:
    python scripts/build_master_catalog.py [--repo <repo_root>]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------
# King mmc4 reader
# ---------------------------------------------------------------------------

def read_king_mmc4(path: Path) -> pd.DataFrame:
    """Read King 2024 mmc4 TF catalog; handle Elsevier header offset."""
    raw = pd.read_excel(path, header=None)
    header_row = None
    for i in range(min(len(raw), 8)):
        vals = [str(x) for x in raw.iloc[i].tolist()[:10]]
        if "Gene ID" in vals and "Human Best Blast Hit" in vals:
            header_row = i
            break
    if header_row is None:
        raise ValueError(f"mmc4 header row not found in {path}")
    df = pd.DataFrame(
        raw.iloc[header_row + 1:].values,
        columns=raw.iloc[header_row].tolist(),
    )
    df = df.dropna(subset=["Gene ID"]).reset_index(drop=True)
    df["Gene ID"] = df["Gene ID"].astype(str).str.strip()
    return df


def _clean(val, max_len: int = 80) -> str:
    s = str(val).strip()
    return "" if s.lower() in ("nan", "none", "") or len(s) > max_len else s


def extract_king_entries(df: pd.DataFrame) -> pd.DataFrame:
    """Pull relevant columns from King mmc4 into a normalized frame."""
    cols = df.columns.tolist()

    # Identify column names flexibly
    def find_col(*keywords) -> str | None:
        for k in keywords:
            for c in cols:
                if k.lower() in str(c).lower():
                    return c
        return None

    id_col = find_col("gene id")
    name_col = find_col("genbank", "planarian gene name", "gene name")
    family_col = find_col("tf family", "family")
    tf_flag_col = find_col("tf?")
    fstf_col = find_col("fstf?")
    blast_col = find_col("human best blast", "blast hit")

    records = []
    for _, r in df.iterrows():
        v6 = _clean(r.get(id_col or "Gene ID", ""))
        if not v6 or not v6.startswith("dd_Smed_v6"):
            continue
        records.append({
            "v6_id": v6,
            "gene_name_king": _clean(r.get(name_col, "")) if name_col else "",
            "tf_family_king": _clean(r.get(family_col, "")) if family_col else "",
            "is_tf_king": str(r.get(tf_flag_col, "")).strip().upper() == "TF"
            if tf_flag_col else False,
            "is_fstf": str(r.get(fstf_col, "")).strip().lower() == "yes"
            if fstf_col else False,
            "human_ortholog_king": _clean(r.get(blast_col, ""), 60)
            if blast_col else "",
        })
    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Perez MOESM5 reader (pre-processed CSV)
# ---------------------------------------------------------------------------

def read_perez_summary(path: Path) -> pd.DataFrame:
    """Load the pre-processed Perez TF summary CSV."""
    df = pd.read_csv(path, dtype=str).fillna("")
    # Standardize column names
    df = df.rename(columns={
        "v6_id": "v6_id",
        "tf_class": "tf_class_perez",
        "tf_family": "tf_family_perez",
        "pfam_domain": "pfam_domain",
        "human_ortholog": "human_ortholog_perez",
    })
    # Only keep rows with a v6 mapping
    df = df[df["v6_id"].str.startswith("dd_Smed_v6")].copy()
    df["v6_id"] = df["v6_id"].str.strip()
    return df


# ---------------------------------------------------------------------------
# Main merge
# ---------------------------------------------------------------------------

def build(repo: Path) -> int:
    king_dir = repo / "datasets" / "raw" / "Supplementary_Data_ King_2024"
    perez_csv = repo / "projects" / "NeuralTF" / "data" / "perez_tf_summary.csv"
    out_path = repo / "projects" / "NeuralTF" / "data" / "master_tf_catalog.csv"

    # Locate mmc4
    mmc4 = king_dir / "1-s2.0-S2211124724001712-mmc4.xlsx"
    if not mmc4.exists():
        for p in sorted(king_dir.iterdir()):
            if p.suffix.lower() == ".xlsx" and "mmc4" in p.stem.lower():
                mmc4 = p
                break
    if not mmc4.exists():
        raise FileNotFoundError(f"King mmc4 not found in {king_dir}")

    # Load King mmc5 for RNAi flag — parses short dd IDs (dd11150) to v6
    mmc5 = king_dir / "1-s2.0-S2211124724001712-mmc5.xlsx"
    rnai_ids: set[str] = set()
    if mmc5.exists():
        raw5 = pd.read_excel(mmc5, header=None)
        # mmc5 col-0 has entries like "dd11150", "dd22163 (UNCX)", "pax2b", etc.
        import re as _re
        _DD_RE = _re.compile(r"dd\w*?(\d{4,6})")
        for _, row in raw5.iterrows():
            val = str(row.iloc[0]).strip()
            if not val or val == "nan":
                continue
            # Collect full v6 IDs matching this short ID
            m = _DD_RE.search(val)
            if m:
                rnai_ids.add(val.split(" ")[0].strip())   # store short form too
                rnai_ids.add(m.group(0).strip())           # e.g. "dd11150"
            elif val.startswith("dd_Smed_v6"):
                rnai_ids.add(val.split(" ")[0])
        print(f"  RNAi short IDs from mmc5: {len(rnai_ids)}")


    print(f"Reading King mmc4: {mmc4.name}")
    king_raw = read_king_mmc4(mmc4)
    king = extract_king_entries(king_raw)
    print(f"  King entries with v6 IDs: {len(king)}")

    print(f"Reading Perez summary: {perez_csv.name}")
    if not perez_csv.exists():
        raise FileNotFoundError(
            f"Perez summary CSV not found: {perez_csv}\n"
            "Run: python projects/NeuralTF/scripts/preprocess_perez.py"
        )
    perez = read_perez_summary(perez_csv)
    print(f"  Perez entries with v6 IDs: {len(perez)}")

    # Outer merge on v6_id
    merged = pd.merge(
        king[["v6_id", "gene_name_king", "tf_family_king", "is_tf_king",
              "is_fstf", "human_ortholog_king"]],
        perez[["v6_id", "tf_class_perez", "tf_family_perez", "pfam_domain",
               "human_ortholog_perez"]].drop_duplicates("v6_id"),
        on="v6_id",
        how="outer",
    ).fillna("")

    # Source flag
    in_king = merged["gene_name_king"] != ""
    in_perez = merged["tf_class_perez"] != ""
    merged["in_king"] = in_king
    merged["in_perez"] = in_perez
    merged["source"] = "perez_only"
    merged.loc[in_king & ~in_perez, "source"] = "king_only"
    merged.loc[in_king & in_perez, "source"] = "both"

    # Unified gene_name: prefer King, fall back to Perez human_ortholog
    merged["gene_name"] = merged["gene_name_king"].where(
        merged["gene_name_king"] != "", merged["human_ortholog_perez"]
    )
    # Unified human_ortholog: prefer Perez, fall back to King
    merged["human_ortholog"] = merged["human_ortholog_perez"].where(
        merged["human_ortholog_perez"] != "", merged["human_ortholog_king"]
    )
    # RNAi flag: rnai_ids has short IDs like 'dd11150'
    # v6_id looks like 'dd_Smed_v6_11150_0_1' — extract numeric portion to match
    import re as _re2
    _NUM_RE = _re2.compile(r"_(\d{4,6})_")
    def _short_from_v6(v6: str) -> str:
        m = _NUM_RE.search(v6)
        return f"dd{m.group(1)}" if m else ""

    merged["has_rnai"] = merged["v6_id"].apply(
        lambda v6: (v6 in rnai_ids)
        or (_short_from_v6(str(v6)) in rnai_ids)
        or any(str(v6).endswith(rid.lstrip("dd")) for rid in rnai_ids)
    )


    # Select and reorder final columns
    out_cols = [
        "v6_id", "gene_name", "tf_family_king", "tf_class_perez",
        "tf_family_perez", "pfam_domain", "human_ortholog",
        "is_fstf", "has_rnai", "in_king", "in_perez", "source",
    ]
    out = merged[out_cols].sort_values("v6_id").reset_index(drop=True)

    # Deduplicate on v6_id (keep 'both' over singles, then first alphabetically)
    source_priority = {"both": 0, "king_only": 1, "perez_only": 2}
    out["_sp"] = out["source"].map(source_priority).fillna(3)
    out = (
        out.sort_values(["v6_id", "_sp"])
        .drop_duplicates(subset="v6_id", keep="first")
        .drop(columns=["_sp"])
        .reset_index(drop=True)
    )

    # Ensure boolean columns are properly typed (outer merge + fillna may leave strings)
    out["is_fstf"] = out["is_fstf"].apply(
        lambda x: bool(x) if isinstance(x, bool) else str(x).lower() == "true"
    )
    out["has_rnai"] = out["has_rnai"].apply(
        lambda x: bool(x) if isinstance(x, bool) else str(x).lower() == "true"
    )
    out["in_king"] = out["in_king"].apply(
        lambda x: bool(x) if isinstance(x, bool) else str(x).lower() == "true"
    )
    out["in_perez"] = out["in_perez"].apply(
        lambda x: bool(x) if isinstance(x, bool) else str(x).lower() == "true"
    )


    # ---------- Assertions ---------------------------------------------------
    king_v6_count = king["v6_id"].nunique()
    assert len(out) >= king_v6_count, (
        f"Master catalog ({len(out)} rows) has fewer entries than King alone "
        f"({king_v6_count}); merge error."
    )
    assert out["v6_id"].nunique() == len(out), "Duplicate v6_ids in master catalog"
    assert "source" in out.columns

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)

    n_both = (out["source"] == "both").sum()
    n_king_only = (out["source"] == "king_only").sum()
    n_perez_only = (out["source"] == "perez_only").sum()
    n_fstf = out["is_fstf"].sum()
    n_rnai = out["has_rnai"].sum()

    print(f"\n[OK] Master TF Catalog written: {out_path}")
    print(f"  Total unique v6 IDs : {len(out):,}")
    print(f"  Both King + Perez   : {n_both:,}")
    print(f"  King only           : {n_king_only:,}")
    print(f"  Perez only          : {n_perez_only:,}")
    print(f"  FSTF flagged        : {n_fstf:,}")
    print(f"  RNAi screen targets : {n_rnai:,}")

    return 0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=ROOT,
                   help="Repository root (default: parent of this file)")
    args = p.parse_args()
    return build(args.repo.resolve())


if __name__ == "__main__":
    sys.exit(main())
