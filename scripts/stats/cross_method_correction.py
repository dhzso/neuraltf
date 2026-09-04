#!/usr/bin/env python
"""Multiple testing correction for cross-method consensus.

Applies Bonferroni and FDR (Benjamini-Hochberg) corrections when
claiming consensus across 3 ranking methods (fixed, centered, uniform).

Usage:
    python scripts/stats/cross_method_correction.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
FIG_DIR = REPO / "projects" / "NeuralTF" / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


def bonferroni_correction(pvalues, alpha=0.05):
    """Bonferroni correction."""
    m = len(pvalues)
    adjusted = np.minimum(np.array(pvalues) * m, 1.0)
    return adjusted


def benjamini_hochberg(pvalues, alpha=0.05):
    """Benjamini-Hochberg FDR correction (step-up, monotone enforced)."""
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    # enforce monotonicity from the largest p downwards
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.zeros(n)
    adjusted[order] = np.minimum(ranked, 1.0)
    return adjusted


def _binom_p(k, n, p):
    """Compute binomial test p-value compatible with SciPy <1.12 and >=1.12."""
    if hasattr(stats, "binomtest"):
        return float(stats.binomtest(k, n, p, alternative="greater").pvalue)
    return float(stats.binom_test(k, n, p, alternative="greater"))


def load_method_top10():
    """Load top-10 from each method."""
    methods = {}

    centered_path = RESULTS_DIR / "dirichlet_centered_full_rank.csv"
    if centered_path.exists():
        df = pd.read_csv(centered_path).drop_duplicates(subset="gene_id", keep="first")
        gene_col = "gene_id" if "gene_id" in df.columns else "gene_id_v6"
        score_col = "composite_score" if "composite_score" in df.columns \
            else "dirichlet_median_score"
        methods["centered"] = set(df.sort_values(score_col, ascending=False)[gene_col].head(10).values)

    uniform_path = RESULTS_DIR / "dirichlet_uniform_full_rank.csv"
    if uniform_path.exists():
        df = pd.read_csv(uniform_path).drop_duplicates(subset="gene_id", keep="first")
        gene_col = "gene_id" if "gene_id" in df.columns else "gene_id_v6"
        score_col = "composite_score" if "composite_score" in df.columns \
            else "uniform_median_score"
        methods["uniform"] = set(df.sort_values(score_col, ascending=False)[gene_col].head(10).values)

    fixed_path = RUN_DIR / "rank.csv"
    if fixed_path.exists():
        df = pd.read_csv(fixed_path)
        gene_col = "gene_id" if "gene_id" in df.columns else "gene_id_v6"
        score_col = "composite_score" if "composite_score" in df.columns \
            else "integrated_score"
        methods["fixed"] = set(
            df.sort_values(score_col, ascending=False)[gene_col].head(10).values
        )

    return methods


def main():
    print("=== Cross-Method Consensus Correction ===")

    methods = load_method_top10()
    print(f"Methods loaded: {list(methods.keys())}")

    if len(methods) < 2:
        print("Error: fewer than 2 methods found")
        return 1

    all_genes = set()
    for s in methods.values():
        all_genes.update(s)
    all_genes = sorted(all_genes)

    n_methods = len(methods)
    n_genes = len(all_genes)

    consensus_data = []
    for gene in all_genes:
        methods_present = [m for m in methods if gene in methods[m]]
        k = len(methods_present)
        p_obs = k / n_methods
        p_random = 1.0 / n_methods

        binom_p = _binom_p(k, n_methods, p_random)
        consensus_data.append({
            "gene_id": gene,
            "n_methods_present": k,
            "methods": ",".join(sorted(methods_present)),
            "p_binom": binom_p,
            "is_consensus": k >= 2,
        })

    df = pd.DataFrame(consensus_data)
    pvals = df["p_binom"].values

    bonf_p = bonferroni_correction(pvals)
    fdr_p = benjamini_hochberg(pvals)

    df["p_bonferroni"] = bonf_p
    df["p_fdr_bh"] = fdr_p
    df["significant_bonferroni"] = bonf_p < 0.05
    df["significant_fdr"] = fdr_p < 0.05

    df = df.sort_values("p_binom")
    df = df.reset_index(drop=True)
    df.index = df.index + 1
    df.index.name = "rank"

    print(f"\nTotal genes in any method top-10: {n_genes}")
    print(f"Genes in >= 2 methods: {(df['n_methods_present'] >= 2).sum()}")
    print(f"Genes in all 3 methods: {(df['n_methods_present'] == 3).sum()}")

    print(f"\nSignificant after Bonferroni (p<0.05): {df['significant_bonferroni'].sum()}")
    print(f"Significant after BH-FDR (p<0.05): {df['significant_fdr'].sum()}")

    print(f"\nTop consensus genes:")
    for _, row in df.head(15).iterrows():
        sig_b = "Y" if row["significant_bonferroni"] else " "
        sig_f = "Y" if row["significant_fdr"] else " "
        print(f"  {row['gene_id']:>30s}  methods={row['n_methods_present']}/{n_methods}  "
              f"p={row['p_binom']:.4e}  Bonf={sig_b}  FDR={sig_f}")

    out_path = RESULTS_DIR / "cross_method_significance.json"
    output = {
        "n_methods": n_methods,
        "method_names": list(methods.keys()),
        "n_genes_any_method": n_genes,
        "n_consensus_2plus": int((df["n_methods_present"] >= 2).sum()),
        "n_consensus_all3": int((df["n_methods_present"] == 3).sum()),
        "n_significant_bonferroni": int(df["significant_bonferroni"].sum()),
        "n_significant_fdr": int(df["significant_fdr"].sum()),
        "genes": df.to_dict(orient="records"),
    }
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved: {out_path}")

    csv_path = RESULTS_DIR / "cross_method_significance.csv"
    df.to_csv(csv_path)
    print(f"Saved: {csv_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
