"""Preprocess Perez 2025 supplementary data for NeuralTF pipeline.

Reads MOESM5 (TF classification + orthologs) and produces a summary CSV
that the pipeline can load directly.

Output: data/perez_tf_summary.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

RAW = ROOT / "datasets" / "raw" / "Supplementary_Data_ Perez_2025"
OUT = ROOT / "projects" / "NeuralTF" / "data" / "perez_tf_summary.csv"


def build() -> int:
    moesm5 = RAW / "41467_2025_65712_MOESM5_ESM.xlsx"
    if not moesm5.exists():
        print(f"  (missing {moesm5}, skipping)")
        return 0

    print(f"Loading {moesm5.name} ...")
    df = pd.read_excel(moesm5, sheet_name=0, dtype=str)
    print(f"  {len(df)} rows, {len(df.columns)} columns")

    cols = df.columns.tolist()
    gene_col = cols[0]  # 'gene' (h1SMcG IDs)

    # Find key columns by pattern matching
    tf_class_col = next((c for c in cols if "TF Class" in c and "Perez" in c), None)
    rbh_col = next((c for c in cols if "1:1" in c and "v6" in c.lower()), None)
    pfam_col = next((c for c in cols if "PFAM domain name" in c and "TF" in c), None)
    family_col = next((c for c in cols if "assigned TF family" in c), None)
    go_col = next((c for c in cols if c.startswith("Gene Ontologies") and "EGGNOG" in c), None)
    ortholog_col = next((c for c in cols if "Collapsed" in c and "Orthofinder" in c), None)
    module_col = next((c for c in cols if "module of TF" in c), None)

    print(f"  TF class column: {tf_class_col}")
    print(f"  RBH v6 column: {rbh_col}")

    records = []
    for _, r in df.iterrows():
        h1 = str(r.get(gene_col, "")).strip()
        if not h1 or h1 == "nan":
            continue

        v6 = str(r.get(rbh_col, "")).strip() if rbh_col else ""
        if v6 == "nan":
            v6 = ""

        tf_class = str(r.get(tf_class_col, "")).strip() if tf_class_col else ""
        if tf_class == "nan":
            tf_class = ""

        pfam = str(r.get(pfam_col, "")).strip() if pfam_col else ""
        if pfam == "nan":
            pfam = ""

        family = str(r.get(family_col, "")).strip() if family_col else ""
        if family == "nan":
            family = ""

        go = str(r.get(go_col, "")).strip() if go_col else ""
        if go == "nan":
            go = ""

        ortholog = str(r.get(ortholog_col, "")).strip() if ortholog_col else ""
        if ortholog == "nan":
            ortholog = ""

        module = str(r.get(module_col, "")).strip() if module_col else ""
        if module == "nan":
            module = ""

        records.append({
            "h1smcg_id": h1,
            "v6_id": v6,
            "tf_class": tf_class,
            "tf_family": family,
            "pfam_domain": pfam,
            "go_terms": go,
            "human_ortholog": ortholog,
            "tf_module": module,
        })

    out_df = pd.DataFrame(records)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT, index=False)

    n_tf = (out_df["tf_class"] != "").sum()
    n_v6 = (out_df["v6_id"] != "").sum()
    print(f"\nSaved {len(out_df)} genes to {OUT}")
    print(f"  With TF class: {n_tf}")
    print(f"  With v6 ID: {n_v6}")
    return 0


if __name__ == "__main__":
    sys.exit(build())
