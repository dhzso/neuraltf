#!/usr/bin/env python
"""Parse King TF Atlas (mmc7.xlsx) -> structured TSV with name->v6 mapping.

3 compartments: G0 Progenitor, X1 (neoblasts), X1 Major Tissue.
Each has 3 aligned sheets: TF Atlas, Pvalues, Log2FC.
Gene names are short aliases -> mapped to v6 IDs via multiple strategies.

Usage: python scripts/build_king_atlas.py
       python -m scripts.build_king_atlas
"""
import re
import pandas as pd
from pathlib import Path

# Repo-relative paths
REPO_ROOT = Path(__file__).resolve().parent.parent
RAW = REPO_ROOT / "datasets" / "raw"
DATA = REPO_ROOT / "projects" / "NeuralTF" / "data"

MMC7 = RAW / "Supplementary_Data_ King_2024" / "1-s2.0-S2211124724001712-mmc7.xlsx"
MMC4 = RAW / "Supplementary_Data_ King_2024" / "1-s2.0-S2211124724001712-mmc4.xlsx"
BRIDGE = DATA / "bridge.csv"

SHEET_NAMES = {
    "G0 Progenitor": ("G0 Progenitor TF Atlas", "G0 Progenitor Pvalues", "G0 Progenitor Log2FC"),
    "X1": ("X1 TF Atlas", "X1 Pvalues", "X1 Log2FC"),
    "X1 Major Tissue": ("X1 Major Tissue Atlas", "X1 Major Tissue Atlas Pvalues", "X1 Major Tissue Atlas Log2F"),
}


def build_name_lookup():
    goods = {}  # lower(name) -> (v6_id, gene_name)

    # From TF catalog GenBank descriptions
    cat = pd.read_excel(MMC4, sheet_name="TF")
    for _, row in cat.iterrows():
        v6 = str(row["Gene ID"]).strip()
        desc = str(row["Planarian GenBank Gene Name"])
        desc_lower = desc.lower()

        m = re.search(r'^([^(]+?)\s*\(not deposited\)', desc)
        if m:
            name = m.group(1).strip()
            goods[name.lower()] = (v6, name)

        m = re.search(r'\(([^()]+)\)\s*(?:mRNA|cds|sequence)', desc, re.IGNORECASE)
        if m:
            name = m.group(1).strip()
            goods[name.lower()] = (v6, name)

        m = re.search(r'(\b[a-zA-Z][a-zA-Z0-9/_-]*)\s+(?:mRNA|cds)', desc, re.IGNORECASE)
        if m:
            name = m.group(1).strip()
            if name.lower() not in goods:
                goods[name.lower()] = (v6, name)

        parts = v6.split("_")
        if len(parts) > 3:
            num = parts[3]
            key = f"dd{num}"
            if key not in goods:
                goods[key] = (v6, key)

    if BRIDGE.exists():
        bridge = pd.read_csv(BRIDGE)
        for _, row in bridge.iterrows():
            name = str(row["gene_name"]).strip()
            v6 = str(row["v6_id"]).strip()
            if name and name != "nan" and v6 and v6 != "nan":
                goods[name.lower()] = (v6, name)
        for _, row in bridge.iterrows():
            for vid_col in ["v6_id", "v4_id"]:
                vid = str(row[vid_col]).strip()
                if vid and vid != "nan" and vid.startswith("dd_Smed_v"):
                    parts = vid.split("_")
                    if len(parts) > 3:
                        num = parts[3]
                        key = f"dd{num}"
                        if key not in goods:
                            goods[key] = (str(row["v6_id"]).strip(), key)

    manual = {
        "gfi1b": "dd_Smed_v6_14824_0_1",
        "runt-1": "dd_Smed_v6_16124_0_1",
        "tcf1": "dd_Smed_v6_13056_0_1",
        "tbx2/3b": "dd_Smed_v6_6470_0_1",
        "sox2": "dd_Smed_v6_8104_0_1",
        "gata4/5/6-3": "dd_Smed_v6_58909_0_1",
        "lhx2/9": "dd_Smed_v6_15144_0_1",
        "scratch": "dd_Smed_v6_18952_0_1",
    }
    for name, v6 in manual.items():
        if name.lower() not in goods:
            goods[name.lower()] = (v6, name)

    for k in [k for k in goods if "/" in k]:
        v = goods[k]
        for sep in ["-", "_"]:
            alt = k.replace("/", sep)
            if alt not in goods:
                goods[alt] = v

    return goods


def _clean_tf_name(raw):
    raw = raw.strip()
    if not raw or raw.startswith("pvalue") or raw.startswith("validated"):
        return []
    m = re.match(r"^(.+?)\s*\(dd(\d+)\)$", raw)
    if m:
        return [(m.group(1).strip(), m.group(2))]
    if re.match(r"^dd\d+$", raw):
        return [(raw, raw[2:])]
    return [(raw, None)]


def _map_to_v6(name, dd_id, goods):
    if dd_id:
        key = f"dd{dd_id}"
        if key in goods:
            return goods[key]
    key = name.lower()
    if key in goods:
        return goods[key]
    key2 = re.sub(r"[-_]\d+$", "", key)
    if key2 != key and key2 in goods:
        return goods[key2]
    key3 = key.replace("-", "_")
    if key3 != key and key3 in goods:
        return goods[key3]
    return (None, None)


def parse_compartment(compartment, goods):
    atlas_sheet, pval_sheet, fc_sheet = SHEET_NAMES[compartment]
    df_atlas = pd.read_excel(MMC7, sheet_name=atlas_sheet, header=None)
    df_pval = pd.read_excel(MMC7, sheet_name=pval_sheet, header=None)
    df_fc = pd.read_excel(MMC7, sheet_name=fc_sheet, header=None)

    is_major = (compartment == "X1 Major Tissue")
    if is_major:
        data_start = 3
        tf_col_start = 1
    else:
        hdr = None
        for i in range(min(10, len(df_atlas))):
            val = str(df_atlas.iloc[i, 0]).strip().lower() if pd.notna(df_atlas.iloc[i, 0]) else ""
            if "subcluster" in val:
                hdr = i
                break
        if hdr is None:
            print(f"  SKIP {compartment}: no header row found")
            return []
        data_start = hdr + 1
        tf_col_start = 3

    records = []
    for idx in range(data_start, len(df_atlas)):
        row = df_atlas.iloc[idx]
        subcluster = str(row.iloc[0]).strip() if pd.notna(row.iloc[0]) else ""
        if not subcluster or subcluster.lower().startswith("table") or subcluster.lower().startswith("log"):
            continue

        if is_major:
            cell_type = ""
            fincher_cl = ""
        else:
            cell_type = str(row.iloc[1]).strip() if pd.notna(row.iloc[1]) else ""
            fincher_cl = str(row.iloc[2]).strip() if pd.notna(row.iloc[2]) else ""

        for c in range(tf_col_start, df_atlas.shape[1]):
            raw = row.iloc[c]
            if pd.isna(raw):
                continue
            raw_str = str(raw).strip()
            if raw_str.startswith("pvalue") or raw_str.startswith("validated"):
                continue

            pval = None
            log2fc = None
            if idx < len(df_pval) and c < df_pval.shape[1]:
                pv = df_pval.iloc[idx, c]
                if pd.notna(pv):
                    pval = float(pv)
            if idx < len(df_fc) and c < df_fc.shape[1]:
                fv = df_fc.iloc[idx, c]
                if pd.notna(fv):
                    log2fc = float(fv)

            for name, dd_id in _clean_tf_name(raw_str):
                v6_id, gene_name = _map_to_v6(name, dd_id, goods)
                if v6_id is None:
                    continue
                if str(pval) == "nan":
                    pval = None
                records.append({
                    "v6_id": v6_id,
                    "gene_name": gene_name or "",
                    "compartment": compartment,
                    "subcluster": subcluster,
                    "cell_type": cell_type,
                    "fincher_cluster": fincher_cl,
                    "log2fc": log2fc,
                    "pval": pval,
                })

    return records


def main():
    print("Building King TF Atlas...")
    goods = build_name_lookup()
    print(f"  Name->v6 lookup: {len(goods)} entries")

    compartmentables = ["G0 Progenitor", "X1", "X1 Major Tissue"]
    all_records = []
    for comp in compartmentables:
        print(f"  Parsing {comp}...")
        recs = parse_compartment(comp, goods)
        print(f"    {len(recs)} mapped records")
        all_records.extend(recs)

    if not all_records:
        print("  No records found!")
        return

    df = pd.DataFrame(all_records)
    df = df.drop_duplicates(subset=["v6_id", "compartment", "subcluster"])

    out_path = DATA / "king_atlas.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  Wrote {len(df)} records to {out_path}")

    for comp in compartmentables:
        sub = df[df["compartment"] == comp]
        print(f"  {comp}: {len(sub)} records, {sub['v6_id'].nunique()} unique TFs")

    print(f"  Neural subclusters: {df['subcluster'].astype(str).str.startswith('neural').sum()} neural rows")


if __name__ == "__main__":
    main()