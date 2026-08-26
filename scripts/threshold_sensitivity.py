#!/usr/bin/env python
"""Threshold sensitivity analysis for NeuralTF pipeline.

Tests stability of top-10 candidates across key threshold parameters.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"


def load_baseline_top10() -> set[str]:
    """Load baseline fixed-weight top-10 gene IDs."""
    path = RESULTS_DIR / "top10_neural_tfs_prioritized.csv"
    if not path.exists():
        raise FileNotFoundError(f"Baseline top10 not found: {path}")
    df = pd.read_csv(path)
    return set(df["gene_id_v6"].tolist())


def run_king_sensitivity(fc_threshold: float) -> set[str]:
    """Re-run King atlas integration with custom log2FC threshold."""
    # This would require re-running pipeline with modified threshold
    # For now, we approximate by filtering existing king_atlas.tsv
    king_path = REPO / "projects" / "NeuralTF" / "data" / "king_atlas.tsv"
    if not king_path.exists():
        return set()

    king = pd.read_csv(king_path, sep="\t")
    neural_mask = king["subcluster"].astype(str).str.startswith("neural")
    neural_df = king[neural_mask & (king["log2fc"] >= fc_threshold)]

    # Get v6_ids of neural-enriched TFs
    return set(neural_df["v6_id"].unique())


def run_cui_sensitivity(neural_mult: float) -> set[str]:
    """Re-run Cui atlas neural enrichment with custom multiplier."""
    cui_path = REPO / "projects" / "NeuralTF" / "data" / "cui_atlas_summary.csv"
    if not cui_path.exists():
        return set()

    cui = pd.read_csv(cui_path)
    # Recompute neural_enriched with custom multiplier
    cui["neural_enriched_custom"] = cui["neural_max"] > (neural_mult * cui["median_expr"])
    return set(cui[cui["neural_enriched_custom"]]["gene_id"].unique())


def compute_jaccard(set1: set, set2: set) -> float:
    """Compute Jaccard index between two sets."""
    if not set1 and not set2:
        return 1.0
    if not set1 or not set2:
        return 0.0
    return len(set1 & set2) / len(set1 | set2)


def main():
    parser = argparse.ArgumentParser(description="NeuralTF threshold sensitivity analysis")
    parser.add_argument("--param", choices=["king_fc", "expr_cap", "cui_neural", "dirichlet_k", "fdr"], required=True)
    parser.add_argument("--values", required=True, help="Comma-separated values to test")
    args = parser.parse_args()

    values = [float(v) for v in args.values.split(",")]
    baseline = load_baseline_top10()

    print(f"=== Threshold Sensitivity: {args.param} ===")
    print(f"Baseline top-10: {len(baseline)} genes")
    print()

    results = []

    if args.param == "king_fc":
        for v in values:
            neural_genes = run_king_sensitivity(v)
            # This is approximate - full pipeline re-run needed for exact top-10
            jac = compute_jaccard(neural_genes, baseline)
            results.append({"threshold": v, "neural_genes": len(neural_genes), "jaccard_vs_baseline": round(jac, 3)})
            print(f"  log2FC={v:.1f}: {len(neural_genes)} neural genes, Jaccard={jac:.3f}")

    elif args.param == "cui_neural":
        for v in values:
            neural_genes = run_cui_sensitivity(v)
            jac = compute_jaccard(neural_genes, baseline)
            results.append({"threshold": v, "neural_genes": len(neural_genes), "jaccard_vs_baseline": round(jac, 3)})
            print(f"  multiplier={v:.1f}: {len(neural_genes)} neural genes, Jaccard={jac:.3f}")

    elif args.param == "expr_cap":
        # Expression cap affects scoring but not candidate set
        # Would need full pipeline re-run
        print("  Expression cap requires full pipeline re-run. Skipping approximation.")
        results = [{"note": "Requires full pipeline re-run"}]

    elif args.param == "dirichlet_k":
        # Dirichlet k affects weighting but not candidate set
        print("  Dirichlet k requires dirichlet_prioritize.py re-run with modified K_DIR.")
        results = [{"note": "Requires Dirichlet script re-run"}]

    elif args.param == "fdr":
        print("  FDR threshold requires pipeline re-run with modified _FDR_THRESHOLD.")
        results = [{"note": "Requires pipeline re-run"}]

    # Save results
    out_path = RESULTS_DIR / f"threshold_sensitivity_{args.param}.json"
    out_path.write_text(json.dumps({
        "parameter": args.param,
        "baseline_top10_count": len(baseline),
        "results": results
    }, indent=2))
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    sys.exit(main())