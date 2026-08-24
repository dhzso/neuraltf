#!/usr/bin/env python
"""Create comprehensive supplementary tables for all methods."""

import pandas as pd
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
IN_DIR = REPO / "projects" / "NeuralTF" / "results"
OUT_DIR = REPO / "projects" / "NeuralTF" / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main() -> int:
    # Load all three method outputs
    fixed = pd.read_csv(IN_DIR / 'top10_neural_tfs_prioritized.csv')
    centered = pd.read_csv(IN_DIR / 'dirichlet_top10_prioritized.csv')
    uniform = pd.read_csv(IN_DIR / 'dirichlet_uniform_top10.csv')

    # S1: Three-method top-10 comparison
    merged = fixed[['gene_id_v6', 'gene_name', 'track', 'composite_score', 'proof_status']].copy()
    merged = merged.rename(columns={'composite_score': 'fixed_composite'})

    centered_sub = centered[['gene_id_v6', 'composite_score', 'dirichlet_median_score']].rename(
        columns={'composite_score': 'centered_composite', 'dirichlet_median_score': 'centered_median'})
    merged = merged.merge(centered_sub, on='gene_id_v6', how='outer')

    uniform_sub = uniform[['gene_id_v6', 'composite_score', 'uniform_median_score']].rename(
        columns={'composite_score': 'uniform_composite', 'uniform_median_score': 'uniform_median'})
    merged = merged.merge(uniform_sub, on='gene_id_v6', how='outer')

    # Add ranks
    for col in ['fixed_composite', 'centered_composite', 'uniform_composite']:
        merged[f'{col.split("_")[0]}_rank'] = merged[col].rank(ascending=False, method='min').astype('Int64')

    # Sort by fixed composite
    merged = merged.sort_values('fixed_composite', ascending=False).reset_index(drop=True)
    merged.to_csv(IN_DIR / 'supplementary_table_S1_method_comparison.csv', index=False)
    print('Created supplementary_table_S1_method_comparison.csv')

    # S2: Full fixed-weight ranked list (all 249)
    rank = pd.read_csv(REPO / 'projects' / 'NeuralTF' / 'runs' / 'pipeline_run' / 'rank.csv')
    rank_neural = pd.read_csv(REPO / 'projects' / 'NeuralTF' / 'runs' / 'pipeline_run' / 'rank_neural.csv')

    # Add composite scores for all 249
    fixed_top10 = pd.read_csv(IN_DIR / 'top10_neural_tfs_prioritized.csv')
    # Use the prioritize script to get composite for all
    from bioforge.projects.neuraltf.prioritize import (
        map_v6_to_v4, prepare_candidates, attach_v4, merge_annotations,
        summarize_annotations, assign_tracks, compute_composite,
    )
    bridge = pd.read_csv(REPO / 'projects' / 'NeuralTF' / 'data' / 'bridge.csv', dtype=str)
    ann_path = REPO / 'datasets' / 'processed' / 'planmine_annotations.parquet'
    ann = pd.read_parquet(ann_path) if ann_path.exists() else pd.DataFrame()

    mapping = map_v6_to_v4(bridge)
    ann_sum = summarize_annotations(ann) if not ann.empty else pd.DataFrame()
    cand = prepare_candidates(rank, mmc4=pd.DataFrame())
    cand = attach_v4(cand, mapping)
    if not ann_sum.empty:
        cand = merge_annotations(cand, ann_sum)
    cand = compute_composite(cand)
    cand.to_csv(IN_DIR / 'supplementary_table_S2_fixed_all249.csv', index=False)
    print('Created supplementary_table_S2_fixed_all249.csv')

    # S3: Full centered Dirichlet ranked list (all 99 candidates)
    centered_full = pd.read_csv(IN_DIR / 'dirichlet_centered_full_rank.csv')
    centered_full.to_csv(IN_DIR / 'supplementary_table_S3_centered_all99.csv', index=False)
    print('Created supplementary_table_S3_centered_all99.csv')

    # S4: Full uniform Dirichlet ranked list
    uniform_full = pd.read_csv(IN_DIR / 'dirichlet_uniform_full_rank.csv')
    uniform_full.to_csv(IN_DIR / 'supplementary_table_S4_uniform_all99.csv', index=False)
    print('Created supplementary_table_S4_uniform_all99.csv')

    # S5: Full 249 uniform Dirichlet
    uniform_249 = pd.read_csv(IN_DIR / 'dirichlet_uniform_all249_full_rank.csv')
    uniform_249.to_csv(IN_DIR / 'supplementary_table_S5_uniform_all249.csv', index=False)
    print('Created supplementary_table_S5_uniform_all249.csv')

    # S6: FSTF rankings
    fstf_19 = pd.read_csv(IN_DIR / 'fstf_ranked_19_neural.csv')
    fstf_43 = pd.read_csv(IN_DIR / 'fstf_ranked_43_all.csv')
    fstf_74 = pd.read_csv(IN_DIR / 'fstf_ranked_74_catalog.csv')
    fstf_19.to_csv(IN_DIR / 'supplementary_table_S6_fstf_19_neural.csv', index=False)
    fstf_43.to_csv(IN_DIR / 'supplementary_table_S7_fstf_43_all.csv', index=False)
    fstf_74.to_csv(IN_DIR / 'supplementary_table_S8_fstf_74_catalog.csv', index=False)
    print('Created FSTF supplementary tables')

    # S7: 99 vs 249 comparison
    comp = pd.read_csv(IN_DIR / 'dirichlet_uniform_all249_summary.txt')
    comp.to_csv(IN_DIR / 'supplementary_table_S9_99vs249.csv', index=False)
    print('Created supplementary_table_S9_99vs249.csv')

    print('\nAll supplementary tables created in', IN_DIR)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())