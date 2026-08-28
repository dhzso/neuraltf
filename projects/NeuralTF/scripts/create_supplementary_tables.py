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
    rank_neural = pd.read_csv(RUN_DIR / "rank_neural.csv")

    centered = pd.read_csv(IN_DIR / "dirichlet_top10_prioritized.csv")
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
        merged[f"{col.split('_')[0]}_rank"] = merged[col].rank(ascending=False, method="min").astype("Int64")

    merged = merged.sort_values("fixed_composite", ascending=False).reset_index(drop=True)
    merged.to_csv(OUT_DIR / "supplementary_table_S1_method_comparison.csv", index=False)
    print("Created supplementary_table_S1_method_comparison.csv")

    rank.to_csv(OUT_DIR / "supplementary_table_S2_fixed_all249.csv", index=False)
    print("Created supplementary_table_S2_fixed_all249.csv")

    centered_full = pd.read_csv(IN_DIR / "dirichlet_centered_full_rank.csv")
    centered_full.to_csv(OUT_DIR / "supplementary_table_S3_centered_all99.csv", index=False)
    print("Created supplementary_table_S3_centered_all99.csv")

    uniform_full = pd.read_csv(IN_DIR / "dirichlet_uniform_full_rank.csv")
    uniform_full.to_csv(OUT_DIR / "supplementary_table_S4_uniform_all99.csv", index=False)
    print("Created supplementary_table_S4_uniform_all99.csv")

    uniform_249 = pd.read_csv(IN_DIR / "dirichlet_uniform_all249_full_rank.csv")
    uniform_249.to_csv(OUT_DIR / "supplementary_table_S5_uniform_all249.csv", index=False)
    print("Created supplementary_table_S5_uniform_all249.csv")

    tf_neural = pd.read_csv(IN_DIR / "tf_ranked_neural_top19.csv")
    tf_all = pd.read_csv(IN_DIR / "tf_ranked_all_top43.csv")
    tf_catalog = pd.read_csv(IN_DIR / "tf_ranked_catalog_top74.csv")
    tf_neural.to_csv(OUT_DIR / "supplementary_table_S6_tf_neural_top19.csv", index=False)
    tf_all.to_csv(OUT_DIR / "supplementary_table_S7_tf_all_top43.csv", index=False)
    tf_catalog.to_csv(OUT_DIR / "supplementary_table_S8_tf_catalog_top74.csv", index=False)
    print("Created TF ranking supplementary tables")

    print("All supplementary tables created in", OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
