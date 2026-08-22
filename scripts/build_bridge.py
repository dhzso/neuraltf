#!/usr/bin/env python
"""Build v4-v6 gene ID bridge CSV from Rosetta Stone + King TF catalog.

The Rosetta Stone file (`smed_20140614.mapping.rosettastone.2020.txt`) has
three tab-separated columns with a header row:
    ref_id  seq_id  transcriptome_id
`ref_id` is the SMED reference gene id; the same ref appears once per
transcriptome with the matching transcript id in `seq_id`. The bridge is
built by pairing the `dd_Smed_v4` and `dd_Smed_v6` transcript rows per ref,
then cross-referencing the King mmc4 TF catalog (whose "Gene ID" column
holds v6 ids) to attach gene names.

Canonical transcript ids look like `dd_Smed_v4_10007_0_1` /
`dd_Smed_v6_10007_0_1`; the rosetta sometimes lists truncated variants of
the same transcript under a ref, so the most specific (full-prefixed,
`_0_1`-terminated) id is preferred.

Writes the bridge CSV (gene_name, v6_id, v4_id) used by the NeuralTF
pipeline (see bioforge.evidence.gene_mapping.BridgeTable).

Usage:
    python scripts/build_bridge.py [--out projects/NeuralTF/data/bridge.csv]
"""

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

V4_PREFIX = "dd_Smed_v4_"
V6_PREFIX = "dd_Smed_v6_"


def _canonical_seq_ids(seq_ids: pd.Series, prefix: str) -> list[str]:
    """Pick the most specific transcript id per ref from a set of variants.

    Within one transcriptome a ref may list e.g. `dd_Smed_v6_10001_0_1`,
    `dd_Smed_v6_10001_0` and `10001_0_1` for the same transcript. Prefer the
    full-prefixed id ending in `_0_1`; fall back to any full-prefixed id,
    then to any id at all. Returns ids sorted for deterministic output.
    """
    ids = [s for s in seq_ids.astype(str).dropna().tolist() if s and s != "nan"]
    if not ids:
        return []
    full = [s for s in ids if s.startswith(prefix)]
    canonical = [s for s in full if s.endswith("_0_1")]
    choice = canonical or full or ids
    return sorted(set(choice))


def build_bridge(rosetta_path, mmc4_path, out_path):
    rosetta = pd.read_csv(rosetta_path, sep="\t", comment="#")
    if not {"ref_id", "seq_id", "transcriptome_id"} <= set(rosetta.columns):
        raise ValueError(
            "Rosetta file must have header columns ref_id, seq_id, "
            f"transcriptome_id; got {list(rosetta.columns)}"
        )

    v4_rows = rosetta[rosetta["transcriptome_id"] == "dd_Smed_v4"]
    v6_rows = rosetta[rosetta["transcriptome_id"] == "dd_Smed_v6"]

    v4_map = {
        ref: _canonical_seq_ids(g["seq_id"], V4_PREFIX)
        for ref, g in v4_rows.groupby("ref_id")
    }
    v6_map = {
        ref: _canonical_seq_ids(g["seq_id"], V6_PREFIX)
        for ref, g in v6_rows.groupby("ref_id")
    }

    records: list[dict[str, str]] = []
    for ref in set(v4_map) & set(v6_map):
        for v4_id in v4_map[ref]:
            for v6_id in v6_map[ref]:
                records.append({"v6_id": v6_id, "v4_id": v4_id})
    bridge = pd.DataFrame(records)
    print(f"  {len(bridge)} ref-paired v4<->v6 rows")
    if bridge.empty:
        raise SystemExit(
            "No ref_id found in both dd_Smed_v4 and dd_Smed_v6 transcriptomes."
        )
    bridge = bridge.drop_duplicates(subset=["v6_id", "v4_id"])

    tf_cat = pd.read_excel(mmc4_path, sheet_name="TF")
    tf_cat["v6_id"] = tf_cat["Gene ID"].astype(str)

    name_cols = [c for c in tf_cat.columns if "GenBank" in c or "Planarian" in c]
    if not name_cols:
        print("Warning: could not identify GenBank name column in mmc4. "
              "Proceeding without names.")
        tf_cat["gene_name"] = ""
    else:
        def _clean_name(x) -> str:
            if x is None or pd.isna(x) or not isinstance(x, str):
                return ""
            x = x.strip()
            return x if 0 < len(x) < 50 else ""

        tf_cat["gene_name"] = tf_cat[name_cols[0]].apply(_clean_name)

    merged = bridge.merge(
        tf_cat[["gene_name", "v6_id"]], on="v6_id", how="left"
    )
    merged = merged[["gene_name", "v6_id", "v4_id"]].sort_values("v6_id")
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