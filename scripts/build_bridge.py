#!/usr/bin/env python
"""Build v4-v6 gene ID bridge CSV from Rosetta Stone + King TF catalog.

Loads the Rosetta Stone mapping file and King mmc4 TF catalog, cross-references
gene names, and writes a bridge CSV (gene_name, v6_id, v4_id) for all bridged
TF genes.

Usage:
    python scripts/build_bridge.py [--out projects/NeuralTF/data/bridge.csv]
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def build_bridge(rosetta_path, mmc4_path, out_path):
    rosetta = pd.read_csv(rosetta_path, sep="\t", comment="#", header=None)
    rosetta.columns = ["v4_id", "v6_id"]
    rosetta["v6_id"] = rosetta["v6_id"].astype(str)
    rosetta["v4_id"] = rosetta["v4_id"].astype(str)

    tf_cat = pd.read_excel(mmc4_path, sheet_name="TF")

    name_cols = [c for c in tf_cat.columns if "GenBank" in c or "Planarian" in c or "Plan" in c]
    if not name_cols:
        print("Warning: could not identify GenBank name column in mmc4. Proceeding without names.")
        tf_cat["gene_name"] = ""
    else:
        tf_cat["gene_name"] = tf_cat[name_cols[0]].astype(str).apply(
            lambda x: x if x and x != "nan" and len(x) < 50 else ""
        )

    tf_cat["v6_id"] = tf_cat["Gene ID"].astype(str)

    merged = tf_cat[["gene_name", "v6_id"]].merge(
        rosetta, on="v6_id", how="outer"
    )
    merged = merged[["gene_name", "v6_id", "v4_id"]]
    merged.to_csv(out_path, index=False)
    bridged = (merged["v4_id"].notna()).sum()
    print(f"Bridge written to {out_path}")
    print(f"  {len(merged)} rows, {bridged} with v4 mapping")


def main():
    parser = argparse.ArgumentParser(description="Build v4-v6 bridge CSV")
    parser.add_argument(
        "--rosetta",
        default=str(ROOT / "datasets" / "raw" / "smed_20140614.mapping.rosettastone.2020.txt"),
    )
    parser.add_argument(
        "--mmc4",
        default=str(ROOT / "datasets" / "raw" / "Supplementary_Data_ King_2024"
                    / "1-s2.0-S2211124724001712-mmc4.xlsx"),
    )
    parser.add_argument(
        "--out",
        default=str(ROOT / "projects" / "NeuralTF" / "data" / "bridge.csv"),
    )
    args = parser.parse_args()
    build_bridge(args.rosetta, args.mmc4, args.out)


if __name__ == "__main__":
    main()