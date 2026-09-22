#!/usr/bin/env python
"""Create supplementary tables for the three prioritization methods.

All tables share the WS2 unification guarantees: one row per gene, the
same candidate universe (rank.csv), and composite scores that carry the
identical bonus mask across methods.

Outputs (projects/NeuralTF/results/):
  supplementary_table_S1_method_comparison.csv  — per-gene scores + ranks
      for fixed / centered / uniform composites, with a shared bonus column
  supplementary_table_S2_fixed_all_candidates.csv
  supplementary_table_S3_centered_all_candidates.csv
  supplementary_table_S4_uniform_all_candidates.csv
"""

import os

import pandas as pd
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
IN_DIR = REPO / "projects" / "NeuralTF" / "results"
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
OUT_DIR = REPO / "projects" / "NeuralTF" / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _atomic_write_csv(df: pd.DataFrame, path: Path, **kwargs) -> None:
    """Write CSV via a .tmp intermediate then atomically replace the target."""
    tmp = path.with_suffix(".csv.tmp")
    df.to_csv(tmp, **kwargs)
    os.replace(str(tmp), str(path))


def _dedup(df: pd.DataFrame, score_col: str) -> pd.DataFrame:
    """Guard: exactly one row per gene (first by descending score)."""
    if "gene_id" in df.columns:
        df = df.sort_values(score_col, ascending=False) \
            .drop_duplicates(subset="gene_id", keep="first")
    return df


def main() -> int:
    fixed_path = IN_DIR / "fixed_full_rank.csv"
    if fixed_path.exists():
        fixed = pd.read_csv(fixed_path)
        fixed = _dedup(fixed, "composite_score")
    else:
        fixed = pd.read_csv(RUN_DIR / "rank.csv")
        fixed = _dedup(fixed, "integrated_score")
        if "composite_score" not in fixed.columns:
            fixed["composite_score"] = fixed["integrated_score"]

    centered = pd.read_csv(IN_DIR / "dirichlet_centered_full_rank.csv")
    centered = _dedup(centered, "composite_score")
    uniform = pd.read_csv(IN_DIR / "dirichlet_uniform_full_rank.csv")
    uniform = _dedup(uniform, "composite_score")

    # S1 — method comparison on the shared all-candidate universe.
    # "fixed_composite" is the fixed method's composite (integrated score
    # + the same bonuses as the other two methods), NOT the raw integrated
    # score — all three columns must be the same quantity class.
    # phenotype_confirmed (2026-09-11 ground-truth fix) is carried through
    # so every supplementary table exposes the corrected label.
    base_cols = ["gene_id", "gene_name", "integrated_score", "proof_status"]
    if "phenotype_confirmed" in fixed.columns:
        base_cols.append("phenotype_confirmed")
    merged = fixed[base_cols].copy()
    merged["fixed_composite"] = fixed["composite_score"]

    centered_sub = centered[
        ["gene_id", "composite_score", "dirichlet_median_score"]
    ].rename(columns={
        "gene_id": "gene_id",
        "composite_score": "centered_composite",
        "dirichlet_median_score": "centered_median",
    })
    merged = merged.merge(centered_sub, on="gene_id", how="left")

    uniform_sub = uniform[
        ["gene_id", "composite_score", "uniform_median_score"]
    ].rename(columns={
        "gene_id": "gene_id",
        "composite_score": "uniform_composite",
        "uniform_median_score": "uniform_median",
    })
    merged = merged.merge(uniform_sub, on="gene_id", how="left")

    for col in ["fixed_composite", "centered_composite", "uniform_composite"]:
        if col in merged.columns:
            key = col.split("_")[0]
            merged[f"{key}_rank"] = merged[col].rank(
                ascending=False, method="min").astype("Int64")

    merged = merged.sort_values("fixed_composite", ascending=False).reset_index(drop=True)
    assert merged["gene_id"].is_unique, "S1 method comparison has duplicate genes"
    _atomic_write_csv(merged, OUT_DIR / "supplementary_table_S1_method_comparison.csv", index=False)
    print("Created supplementary_table_S1_method_comparison.csv "
          f"({len(merged)} genes)")

    _atomic_write_csv(fixed, OUT_DIR / "supplementary_table_S2_fixed_all_candidates.csv", index=False)
    print("Created supplementary_table_S2_fixed_all_candidates.csv")

    _atomic_write_csv(centered, OUT_DIR / "supplementary_table_S3_centered_all_candidates.csv", index=False)
    print("Created supplementary_table_S3_centered_all_candidates.csv")

    _atomic_write_csv(uniform, OUT_DIR / "supplementary_table_S4_uniform_all_candidates.csv", index=False)
    print("Created supplementary_table_S4_uniform_all_candidates.csv")

    # Optional TF tables if generated
    for src, dst in [
        ("tf_ranked_neural_top19.csv", "supplementary_table_S5_tf_neural.csv"),
        ("tf_ranked_all_top43.csv", "supplementary_table_S6_tf_all.csv"),
        ("tf_ranked_catalog_top74.csv", "supplementary_table_S7_tf_catalog.csv"),
    ]:
        p = IN_DIR / src
        if p.exists():
            _atomic_write_csv(pd.read_csv(p), OUT_DIR / dst, index=False)
            print(f"Created {dst}")

    print("All supplementary tables created in", OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
