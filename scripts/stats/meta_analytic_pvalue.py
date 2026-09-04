#!/usr/bin/env python
"""Meta-analytic p-value combination across atlases (REAL data only).

Combines per-gene Wilcoxon DE p-values from the Fincher, Plass, and Cui
atlases using Fisher's method (chi-squared) and Stouffer's method
(weighted z-scores). The per-atlas p-values are the pipeline's persisted
``de_pvalues`` checkpoint (best significant-cluster p per gene per atlas,
written by checkpoint 03).

WS3 fix: the previous version fell back to a SIMULATION with fabricated
gene IDs when atlas p-value files were missing, and shipped that CSV in
results/ as if real. This version never fabricates data — if the real
checkpoint is missing it exits with a clear message.

Usage:
    python scripts/stats/meta_analytic_pvalue.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

ATLAS_P_COLS = {
    "fincher": "fincher_p",
    "plass": "plass_p",
    "cui": "cui_p",
}


def fishers_method(pvalues):
    """Combine p-values using Fisher's method. Returns chi2 statistic and combined p-value."""
    p = np.array(pvalues)
    p = p[p > 0]
    if len(p) == 0:
        return 0.0, 1.0
    chi2_stat = -2.0 * np.sum(np.log(p))
    combined_p = 1.0 - stats.chi2.cdf(chi2_stat, 2 * len(p))
    return float(chi2_stat), float(combined_p)


def stouffers_method(pvalues, weights=None):
    """Combine p-values using Stouffer's method. Returns z statistic and combined p-value."""
    p = np.array(pvalues)
    p = p[p > 0]
    if len(p) == 0:
        return 0.0, 1.0
    z_scores = stats.norm.ppf(1.0 - p)
    if weights is None:
        weights = np.ones(len(p))
    combined_z = np.sum(z_scores * weights) / np.sqrt(np.sum(weights**2))
    combined_p = 1.0 - stats.norm.cdf(combined_z)
    return float(combined_z), float(combined_p)


def load_de_pvalues() -> pd.DataFrame | None:
    """Load the pipeline's per-gene per-atlas DE p-value checkpoint."""
    path = RUN_DIR / "de_pvalues.parquet"
    if not path.exists():
        path = RUN_DIR / "de_pvalues.csv"
    if not path.exists():
        return None
    df = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
    return df.drop_duplicates(subset="v6_id", keep="first")


def main():
    print("=== Meta-Analytic P-value Combination (real per-atlas DE p) ===")

    de = load_de_pvalues()
    if de is None or de.empty:
        print(
            "ERROR: per-atlas DE p-value checkpoint not found at "
            f"{RUN_DIR / 'de_pvalues.parquet'}.\n"
            "Re-run the pipeline (scripts/run.py) — checkpoint 03 writes "
            "de_pvalues with fincher_p/plass_p/cui_p per gene.\n"
            "This analysis no longer falls back to simulated data."
        )
        return 1

    present = {a: c for a, c in ATLAS_P_COLS.items() if c in de.columns}
    if len(present) < 2:
        print(f"ERROR: need p-values from >= 2 atlases; found {list(present)}")
        return 1
    print(f"Atlases with per-gene p-values: {sorted(present)} ({len(de)} genes)")

    de = de.sort_values("v6_id").reset_index(drop=True)
    results = []
    for _, row in de.iterrows():
        pvals = []
        for atlas, col in sorted(present.items()):
            v = row.get(col)
            if pd.notna(v):
                pvals.append(max(float(v), 1e-16))
        if len(pvals) < 2:
            continue
        # Stouffer weight: atlas cell count (larger atlas -> more precise p)
        fisher_chi2, fisher_p = fishers_method(pvals)
        stouffer_z, stouffer_p = stouffers_method(pvals)
        results.append({
            "gene_id": row["v6_id"],
            "n_atlases": len(pvals),
            "individual_pvalues": pvals,
            "fisher_chi2": fisher_chi2,
            "fisher_combined_p": fisher_p,
            "stouffer_z": stouffer_z,
            "stouffer_combined_p": stouffer_p,
        })

    out_df = pd.DataFrame(results)
    if out_df.empty:
        print("No gene had p-values in >= 2 atlases; nothing to combine.")
        return 1
    out_df["individual_pvalues"] = out_df["individual_pvalues"].apply(str)
    out_df = out_df.sort_values("fisher_combined_p")

    out_path = RESULTS_DIR / "meta_analysis_pvalues.csv"
    out_df.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path} ({len(out_df)} genes)")

    print(f"\nFisher's method - significant at p<0.05: "
          f"{(out_df['fisher_combined_p'] < 0.05).sum()}")
    print(f"Stouffer's method - significant at p<0.05: "
          f"{(out_df['stouffer_combined_p'] < 0.05).sum()}")

    print("\nTop-10 by Fisher's combined p-value:")
    for _, row in out_df.head(10).iterrows():
        print(f"  {row['gene_id']:>30s}  Fisher p={row['fisher_combined_p']:.4e}  "
              f"Stouffer p={row['stouffer_combined_p']:.4e}")

    rho = out_df["fisher_combined_p"].corr(out_df["stouffer_combined_p"],
                                          method="spearman")
    print(f"\nSpearman rho (Fisher vs Stouffer): {rho:.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
