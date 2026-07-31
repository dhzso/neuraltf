#!/usr/bin/env python
"""Build king_atlas.tsv from King mmc7.xlsx.

Extracts the G0 Progenitor TF Atlas log2FC and pvalue matrices and
writes a tidy TSV. Uses the prebuilt header mappings from the Atlas
sheet for subcluster / cell_type / fincher_cluster annotations.

Usage:
    python scripts/build_king_atlas.py [--mmc7 <path>] [--out <path>]
"""

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def _safe_col(df, col_idx, row_idx):
    try:
        val = df.iloc[row_idx, col_idx]
        return val if not pd.isna(val) else ""
    except IndexError:
        return ""


def build_king_atlas(mmc7_path, out_path):
    xl = pd.ExcelFile(mmc7_path)

    atlas = pd.read_excel(xl, sheet_name="G0 Progenitor TF Atlas", header=4)
    fc = pd.read_excel(xl, sheet_name="G0 Progenitor Log2FC", header=4)
    pv = pd.read_excel(xl, sheet_name="G0 Progenitor Pvalues", header=4)

    n_genes = fc.shape[0]
    n_subs = fc.shape[1] - 1

    rows = []
    for i in range(n_genes):
        gene = str(fc.iloc[i, 0]).strip()
        if not gene or gene == "nan":
            continue
        for j in range(n_subs):
            l2fc_val = fc.iloc[i, j + 1]
            if pd.isna(l2fc_val):
                continue
            try:
                l2fc = float(l2fc_val)
            except (ValueError, TypeError):
                continue

            sub_name = str(atlas.columns[j + 1]) if j + 1 < len(atlas.columns) else ""
            cell_type = str(atlas.iloc[0, j + 1]) if atlas.shape[0] > 0 and j + 1 < atlas.shape[1] else ""
            fincher_cl = str(atlas.iloc[1, j + 1]) if atlas.shape[0] > 1 and j + 1 < atlas.shape[1] else ""

            if cell_type == "nan":
                cell_type = ""
            if fincher_cl == "nan":
                fincher_cl = ""

            pval = None
            try:
                pv_val = pv.iloc[i, j + 1]
                if not pd.isna(pv_val):
                    pval = float(pv_val)
            except (IndexError, ValueError, TypeError):
                pass

            rows.append({
                "v6_id": gene if gene.startswith("dd_Smed_v6_") else "",
                "gene_name": gene if not gene.startswith("dd_Smed_v6_") else "",
                "compartment": "G0 Progenitor",
                "subcluster": sub_name,
                "cell_type": cell_type,
                "fincher_cluster": fincher_cl,
                "log2fc": l2fc,
                "pval": pval,
            })

    df = pd.DataFrame(rows)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"King atlas written to {out_path}")
    print(f"  {len(df)} rows, {df['subcluster'].nunique()} subclusters")


def main():
    parser = argparse.ArgumentParser(description="Build King atlas TSV")
    parser.add_argument("--mmc7")
    parser.add_argument("--out")
    args = parser.parse_known_args()[0]
    mmc7 = args.mmc7 or str(
        ROOT / "datasets" / "raw" / "Supplementary_Data_ King_2024"
        / "1-s2.0-S2211124724001712-mmc7.xlsx")
    out = args.out or str(ROOT / "projects" / "NeuralTF" / "data" / "king_atlas.tsv")
    build_king_atlas(mmc7, out)


if __name__ == "__main__":
    main()