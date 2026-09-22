#!/usr/bin/env python
"""Parse King TF Atlas (mmc7.xlsx) -> structured TSV with name->v6 mapping.

3 compartments: G0 Progenitor, X1 (neoblasts), X1 Major Tissue.
Each has 3 aligned sheets: TF Atlas, Pvalues, Log2FC.
Gene names are short aliases -> mapped to v6 IDs via multiple strategies.

2026-09-13 CORRECTNESS FIXES (senior-review audit):
  1. ROW ALIGNMENT: value sheets are now keyed by SUBCLUSTER NAME, not by
     positional row index. The old positional join was off by one for G0
     and X1 (the TF Atlas sheets carry an extra legend row that the
     Pvalues/Log2FC sheets lack), which attached the NEXT subcluster's
     statistics to every emitted row (~836/1,046 rows carried another
     gene's p-value/log2FC from a different subcluster).
  2. X1 COLUMN OFFSET: the X1 TF Atlas has no Fincher-cluster column —
     its TFs start one column earlier than G0. The old code dropped each
     X1 subcluster's lead TF (e.g. pax6A never entered via X1) and
     stuffed TF names into fincher_cluster. X1 fincher_cluster is now "".
  3. FONT-COLOUR SEMANTICS corrected per the sheets' OWN legend:
       red   (FF0000) = highly significant p (<1E-30 G0 / <1E-100 X1)
       green (00B050) = LESS-significant p (>1E-10 G0 / >1E-50 X1)
       blue          = previously published FSTF — NOT PRESENT in the
                       distributed file (verified: only red/green/theme-1
                       fonts exist; no blue font anywhere).
     The old code mislabelled green as "validated/previously reported
     FSTF". ``previously_reported_fstf`` is therefore DEPRECATED and
     always False; ``sig_color`` keeps red/green with the correct
     meaning (significance tier). A renamed ``less_significant`` column
     is emitted for clarity.
  4. NAME RESOLUTION for phenotype-confirmed ground-truth genes that
     mmc4/bridge cannot resolve (IRX1/IRX2/IRX6/Tbx2/3c/Post-2b/RREB2/
     GCM2 etc. have no GenBank description and no PlanMine official
     name). The curated paper-derived name->v6 map (the same one the
     ground-truth module documents) is now applied BEFORE dropping
     unresolved cells; ``IRX6 (dd18125)``-style disambiguating labels
     are honoured, and bare-number parentheticals ("IRX2 (11500)") are
     parsed.

Usage: python scripts/build_king_atlas.py
       python -m scripts.build_king_atlas
"""
import re
import pandas as pd
from pathlib import Path

import openpyxl

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

# Font colours in the TF Atlas sheets (verified against the distributed
# file AND the sheets' own legend text):
RED_RGB = "FFFF0000"     # highly significant p (< 1E-30 G0 / < 1E-100 X1)
GREEN_RGB = "FF00B050"   # LESS-significant p (> 1E-10 G0 / > 1E-50 X1)
# NOTE: the legend describes BLUE as "previously published FSTF", but no
# blue font exists anywhere in the distributed file (verified cell by
# cell: only red, green, and theme-1 black). previously_reported_fstf is
# therefore not recoverable from this file and stays False.

# Curated paper-derived symbol -> v6 map for genes that neither mmc4's
# GenBank descriptions, the bridge, nor PlanMine can resolve. Every entry
# below is independently verifiable from the repo's own sources:
#   - mmc4 GenBank descriptions ('T-box 2/3c protein', 'fer3l-1 protein',
#     'ascl-2 protein', 'Pax6A mRNA', 'nkx2.2 (not deposited)')
#   - mmc5's own target labels ('post2b (dd9061)')
#   - mmc7's own disambiguating labels ('IRX2 (11500)')
#   - PlanMine official names ('tbx2/3b' -> dd6470, 'glial cells missing'
#     description family for GCM2)
#   - the King 2024 paper figure labels as curated in
#     bioforge.evidence.groundtruth (IRX1/IRX6/RREB2 phenotype ties)
PAPER_SYMBOL_MAP_VERIFIED = {
    "tbx2/3c": "dd_Smed_v6_11693_0_1",   # mmc4: 'T-box 2/3c protein mRNA'
    "fer3l-1": "dd_Smed_v6_8096_0_1",    # mmc4: 'fer3l-1 protein (fer3l-1) mRNA'
    "ascl-2": "dd_Smed_v6_14753_0_1",    # mmc4: 'ascl-2 protein (ascl-2) mRNA'
    "pax6a": "dd_Smed_v6_17726_0_1",     # mmc4: 'clone Smed_03546_V2 Pax6A mRNA'
    "nkx2.2": "dd_Smed_v6_2716_0_1",     # mmc4: 'nkx2.2 (not deposited)'
    "nkx2-like": "dd_Smed_v6_11198_0_1",  # mmc4: 'nkx6-like mRNA' (label alias)
    "post-2b": "dd_Smed_v6_9061_0_1",    # mmc5's own label 'post2b (dd9061)'
    "gcm2": "dd_Smed_v6_7752_0_1",       # paper GCM2 phenotype + PlanMine 'glial cells missing' family
    "irx2": "dd_Smed_v6_11500_0_1",      # mmc7's own label 'IRX2 (11500)'
    "rreb2": "dd_Smed_v6_10103_0_1",     # paper Fig S4D phenotype label (curated set)
    "irx1": "dd_Smed_v6_14656_0_1",      # paper Fig S4F phenotype label (curated set)
    "irx6": "dd_Smed_v6_17352_0_1",      # paper Fig 4E/S8A phenotype label (curated set)
    "p53": "dd_Smed_v6_5563_0_1",        # mmc4: 'SMED-P53 mRNA' (case-insensitive)
    "otxb": "dd_Smed_v6_15516_0_1",      # mmc4: 'otxb-like mRNA'
    "nk4": "dd_Smed_v6_9259_0_1",        # mmc4: 'NK homeobox 4 (nk4) mRNA'
}
# Deliberately NOT mapped (no independent verification available from any
# source in this repo — mmc4 has no description, PlanMine has no name,
# and the paper label carries no dd#### disambiguation):
#   mecom, nkx2-4, ascl4, gata3, hoxb7, prdm1-1, cux-1, dlx, egr-3,
#   atoh8-2, pax2b, fer3l-2 — these cells are skipped (and logged).


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
        # 2026-09-13: 'Tbx2/3b' remapped to dd17143 — the paper's own RNAi
        # row (mmc5: dd17143, marker dd210/GLIPR1 = the Fig 4E/S8G Tbx2/3b
        # phenotype) and mmc4's 'T-box 2/3b protein'. dd6470 (PlanMine
        # alias 'tbx2/3b') is in neither table; see groundtruth.py.
        "tbx2/3b": "dd_Smed_v6_17143_0_1",
        "sox2": "dd_Smed_v6_8104_0_1",
        "gata4/5/6-3": "dd_Smed_v6_58909_0_1",
        "lhx2/9": "dd_Smed_v6_15144_0_1",
        "scratch": "dd_Smed_v6_18952_0_1",
    }
    for name, v6 in manual.items():
        if name.lower() not in goods:
            goods[name.lower()] = (v6, name)

    # Paper-derived verified symbol map (fills mmc4/bridge gaps; loaded
    # LAST so verified sources win where they exist — entries here are
    # only those that could NOT be resolved any other way, or whose
    # paper label differs from the GenBank description).
    for name, v6 in PAPER_SYMBOL_MAP_VERIFIED.items():
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
    """Yield (name, dd_id) pairs from a raw TF-cell string.

    Handles the label dialects found in mmc7:
        'pou4l-1'                       -> symbol only
        'dd22331'                       -> dd token only
        'dd22163 (UNCX)'                -> dd token + symbol
        'SOX2 (dd8104)'                 -> symbol + dd token
        'IRX6 (dd18125)'                -> symbol + EXPLICIT v6 number
        'IRX2 (11500)'                  -> symbol + bare number (no 'dd')
    """
    raw = raw.strip()
    if not raw or raw.lower().startswith("pvalue") or raw.lower().startswith("validated"):
        return []
    # symbol + explicit dd-prefixed number: "SOX2 (dd8104)"
    m = re.match(r"^(.+?)\s*\(dd(\d+)\)$", raw, re.IGNORECASE)
    if m:
        return [(m.group(1).strip(), m.group(2))]
    # symbol + bare number: "IRX2 (11500)" — the number IS the v6 numeric
    m = re.match(r"^(.+?)\s*\((\d+)\)$", raw)
    if m and not re.match(r"^dd", m.group(1), re.IGNORECASE):
        return [(m.group(1).strip(), m.group(2))]
    if re.match(r"^dd\d+$", raw, re.IGNORECASE):
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


def _load_color_lookup():
    """Build {(row, col): color_rgb} for the TF Atlas sheets via openpyxl.

    Returns a dict of {sheet_name: {(row_idx, col_idx): rgb_string}} where
    row/col are 1-based openpyxl coordinates. The distributed mmc7 xlsx
    preserves the red/green font encoding (unlike mmc5).

    2026-09-13: semantics corrected per the sheet legend —
      red   = highly significant p-value
      green = LESS-significant p-value (>1E-10 / >1E-50)
      blue  = previously published FSTF — absent from the distributed file.
    """
    wb = openpyxl.load_workbook(MMC7, data_only=True)
    out = {}
    for sheets in SHEET_NAMES.values():
        sheet = sheets[0]  # the TF Atlas sheet carries the colours
        if sheet not in wb.sheetnames:
            continue
        ws = wb[sheet]
        lookup = {}
        for row in ws.iter_rows(min_row=1, max_row=ws.max_row,
                                max_col=ws.max_column):
            for cell in row:
                if cell.value is None:
                    continue
                c = cell.font.color if cell.font else None
                rgb = c.rgb if (c is not None and c.type == "rgb") else None
                if rgb in (RED_RGB, GREEN_RGB):
                    lookup[(cell.row, cell.column)] = rgb
        out[sheet] = lookup
    return out


def _sheet_subcluster_map(df, sub_col=0):
    """Map subcluster-name -> row index for a sheet (first occurrence wins)."""
    out = {}
    for i in range(len(df)):
        v = df.iloc[i, sub_col]
        if pd.notna(v):
            key = str(v).strip()
            if key and key not in out:
                out[key] = i
    return out


def parse_compartment(compartment, goods, color_lookup=None):
    atlas_sheet, pval_sheet, fc_sheet = SHEET_NAMES[compartment]
    df_atlas = pd.read_excel(MMC7, sheet_name=atlas_sheet, header=None)
    df_pval = pd.read_excel(MMC7, sheet_name=pval_sheet, header=None)
    df_fc = pd.read_excel(MMC7, sheet_name=fc_sheet, header=None)

    # openpyxl colours are 1-based; pandas positions are 0-based
    sheet_colors = (color_lookup or {}).get(atlas_sheet, {})

    # Value-sheet row lookup BY SUBCLUSTER NAME (fixes the old positional
    # off-by-one: the Pvalues/Log2FC sheets lack the atlas sheets' extra
    # legend row, so positional indices never aligned for G0/X1).
    pval_row = _sheet_subcluster_map(df_pval)
    fc_row = _sheet_subcluster_map(df_fc)

    is_major = (compartment == "X1 Major Tissue")
    is_x1 = (compartment == "X1")
    if is_major:
        data_start = 3
        tf_col_start = 1
        # Major Tissue sheets have the SAME shape for atlas and value
        # sheets (verified: 3 legend rows + data from row 3, TFs from
        # col 1 in both) — but we key by subcluster name anyway.
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
        # G0 has a Fincher-clusters column (col 2); X1 does NOT — its TFs
        # start at col 2 directly (verified: X1 row 5 = [0_0, Neural,
        # pax6A, runt-1, ...]). The old code used 3 for both, dropping
        # each X1 subcluster's lead TF and mislabelling TF names as
        # Fincher clusters.
        tf_col_start = 2 if is_x1 else 3

    # Column offset between the atlas sheet and the value sheets: NONE —
    # verified empirically for all three compartments. The X1 Pvalues/
    # Log2FC header rows contain a stale 'Fincher Clusters' label in the
    # col-2 header cell (copied from G0), but the DATA in both sheets
    # starts at the same column as the atlas sheet's TF columns:
    #   G0 atlas:  [subcluster, cell_type, fincher_cluster, TFs at 3]
    #   G0 pv/fc:  [subcluster, cell_type, fincher_cluster, values at 3]
    #   X1 atlas:  [subcluster, cell_type, TFs at 2]
    #   X1 pv/fc:  [subcluster, cell_type, values at 2] (header label lies)
    #   Major:     [subcluster, TFs at 1] (both atlas and values)
    val_col_offset = 0

    records = []
    unresolved = []
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
            # X1 has no Fincher-cluster column; col 2 is a TF (handled
            # by tf_col_start) so fincher_cluster stays empty for X1.
            fincher_cl = "" if is_x1 else (str(row.iloc[2]).strip() if pd.notna(row.iloc[2]) else "")

        # Value-sheet rows keyed by subcluster NAME (not position).
        pv_i = pval_row.get(subcluster)
        fc_i = fc_row.get(subcluster)

        for c in range(tf_col_start, df_atlas.shape[1]):
            raw = row.iloc[c]
            if pd.isna(raw):
                continue
            raw_str = str(raw).strip()
            if raw_str.lower().startswith("pvalue") or raw_str.lower().startswith("validated"):
                continue

            pval = None
            log2fc = None
            vc = c + val_col_offset
            if pv_i is not None and vc < df_pval.shape[1]:
                pv = df_pval.iloc[pv_i, vc]
                if pd.notna(pv):
                    try:
                        pval = float(pv)
                    except (TypeError, ValueError):
                        pval = None
            if fc_i is not None and vc < df_fc.shape[1]:
                fv = df_fc.iloc[fc_i, vc]
                if pd.notna(fv):
                    try:
                        log2fc = float(fv)
                    except (TypeError, ValueError):
                        log2fc = None

            # Font-colour annotations (2026-09-13 corrected semantics):
            #   red   = highly significant p (<1E-30 G0 / <1E-100 X1)
            #   green = LESS-significant p (>1E-10 G0 / >1E-50 X1)
            # previously_reported_fstf is NOT recoverable from this file
            # (the legend's blue font does not exist in the distributed
            # copy) and stays False.
            rgb = sheet_colors.get((idx + 1, c + 1))
            sig_color = "red" if rgb == RED_RGB else ("green" if rgb == GREEN_RGB else "")
            less_significant = sig_color == "green"

            for name, dd_id in _clean_tf_name(raw_str):
                v6_id, gene_name = _map_to_v6(name, dd_id, goods)
                if v6_id is None:
                    unresolved.append((name, subcluster))
                    continue
                records.append({
                    "v6_id": v6_id,
                    "gene_name": gene_name or "",
                    "compartment": compartment,
                    "subcluster": subcluster,
                    "cell_type": cell_type,
                    "fincher_cluster": fincher_cl,
                    "log2fc": log2fc,
                    "pval": pval,
                    "sig_color": sig_color,
                    "less_significant": less_significant,
                    "previously_reported_fstf": False,
                })

    if unresolved:
        from collections import Counter
        cnt = Counter(u[0] for u in unresolved)
        print(f"    [unresolved symbols] {len(unresolved)} cells across "
              f"{len(cnt)} symbols: {dict(cnt.most_common(12))}")

    return records


def main():
    print("Building King TF Atlas...")
    goods = build_name_lookup()
    print(f"  Name->v6 lookup: {len(goods)} entries")

    print("  Reading mmc7 font-colour annotations "
          "(red=highly-significant p, green=less-significant p)...")
    color_lookup = _load_color_lookup()
    n_color_cells = sum(len(v) for v in color_lookup.values())
    print(f"    {n_color_cells} coloured TF cells found across TF Atlas sheets")

    compartmentables = ["G0 Progenitor", "X1", "X1 Major Tissue"]
    all_records = []
    for comp in compartmentables:
        print(f"  Parsing {comp}...")
        recs = parse_compartment(comp, goods, color_lookup=color_lookup)
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

    neural_n = df["subcluster"].astype(str).str.startswith("neural").sum()
    x1_neural_n = (
        (df["compartment"] == "X1") &
        (df["cell_type"].astype(str).str.strip().str.lower() == "neural")
    ).sum()
    print(f"  Neural rows (subcluster prefix 'neural'): {neural_n}")
    print(f"  X1 Neural-cell-type rows: {x1_neural_n}")

    # Sanity check: phenotype-confirmed ground-truth coverage
    import sys
    sys.path.insert(0, str(REPO_ROOT / "src"))
    from bioforge.evidence.groundtruth import PHENOTYPE_CONFIRMED_V6
    covered = PHENOTYPE_CONFIRMED_V6 & set(df["v6_id"])
    print(f"  Phenotype-confirmed genes in atlas: {len(covered)}/{len(PHENOTYPE_CONFIRMED_V6)}")


if __name__ == "__main__":
    main()
