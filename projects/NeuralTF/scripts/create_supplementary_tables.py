#!/usr/bin/env python
"""Create comprehensive supplementary tables for all methods."""

import pandas as pd
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
IN_DIR = REPO / "projects" / "NeuralTF" / "results"
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
OUT_DIR = REPO / "projects" / "NeuralTF" / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main() -> int:
    rank = pd.read_csv(RUN_DIR / "rank.csv")

    centered_path = IN_DIR / "dirichlet_centered_top10.csv"
    if not centered_path.exists():
        centered_path = IN_DIR / "dirichlet_top10_prioritized.csv"
    centered = pd.read_csv(centered_path)

    uniform = pd.read_csv(IN_DIR / "dirichlet_uniform_top10.csv")

    merged = rank[["gene_id", "gene_name", "integrated_score", "proof_status"]].copy()
    merged = merged.rename(columns={"integrated_score": "fixed_composite"})

    centered_sub = centered[["gene_id_v6", "composite_score", "dirichlet_median_score"]].rename(
        columns={"gene_id_v6": "gene_id", "composite_score": "centered_composite",
                 "dirichlet_median_score": "centered_median"})
    merged = merged.merge(centered_sub, on="gene_id", how="outer")

    uniform_sub = uniform[["gene_id_v6", "composite_score", "uniform_median_score"]].rename(
        columns={"gene_id_v6": "gene_id", "composite_score": "uniform_composite",
                 "uniform_median_score": "uniform_median"})
    merged = merged.merge(uniform_sub, on="gene_id", how="outer")

    for col in ["fixed_composite", "centered_composite", "uniform_composite"]:
        if col in merged.columns:
            merged[f"{col.split('_')[0]}_rank"] = merged[col].rank(ascending=False, method="min").astype("Int64")

    merged = merged.sort_values("fixed_composite", ascending=False).reset_index(drop=True)
    merged.to_csv(OUT_DIR / "supplementary_table_S1_method_comparison.csv", index=False)
    print("Created supplementary_table_S1_method_comparison.csv")

    rank.to_csv(OUT_DIR / "supplementary_table_S2_fixed_all_candidates.csv", index=False)
    print("Created supplementary_table_S2_fixed_all_candidates.csv")

    centered_full = pd.read_csv(IN_DIR / "dirichlet_centered_full_rank.csv")
    centered_full.to_csv(OUT_DIR / "supplementary_table_S3_centered_all_candidates.csv", index=False)
    print("Created supplementary_table_S3_centered_all_candidates.csv")

    uniform_full = pd.read_csv(IN_DIR / "dirichlet_uniform_full_rank.csv")
    uniform_full.to_csv(OUT_DIR / "supplementary_table_S4_uniform_all_candidates.csv", index=False)
    print("Created supplementary_table_S4_uniform_all_candidates.csv")

    # Optional TF tables if generated
    for src, dst in [
        ("tf_ranked_neural_top19.csv", "supplementary_table_S5_tf_neural.csv"),
        ("tf_ranked_all_top43.csv", "supplementary_table_S6_tf_all.csv"),
        ("tf_ranked_catalog_top74.csv", "supplementary_table_S7_tf_catalog.csv"),
    ]:
        p = IN_DIR / src
        if p.exists():
            pd.read_csv(p).to_csv(OUT_DIR / dst, index=False)

    print("All supplementary tables created in", OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
