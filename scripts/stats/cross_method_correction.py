#!/usr/bin/env python
"""Cross-method consensus analysis with a VALID randomization null.

2026-09-04 statistical redesign (audit findings M1):

1. The "fixed" arm now loads the ACTUAL published fixed-method shortlist
   (results/top10_neural_tfs_prioritized.csv - composite score + bonus
   mask + dual-track 5+5 selection), not a raw integrated_score top-10
   from rank.csv (the old arm shared only 4/10 genes with the real
   shortlist).
2. The consensus null is the randomization probability that a gene lands
   in a given method's top-10 by chance, p0 = n_top / N_universe (with
   N = 11,672, p0 = 10/11,672 = 8.6e-4) - NOT the previous binomial
   p=1/3 "fair coin" null, which was invalid by ~390x and produced a
   structurally zero-power test (it could never find significance and
   concluded "no consensus" when 3/3 overlap is in fact overwhelming
   evidence, p ~ 1e-9). This now matches the null used by
   overlap_significance.py (10/N).
3. Significance is reported both per-gene (binomial k-of-n_methods at
   p0) and as the global top-10 set overlap (hypergeometric), with
   Bonferroni/BH across the tested gene family.
4. A documented caveat: the three methods share the candidate matrix and
   (by design) the bonus layer, so membership events are positively
   correlated; the binomial p-values are therefore descriptive of
   agreement-strength, and the hypergeometric set-level test is the
   primary consensus statistic.

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

N_TOP = 10  # each method's shortlist size


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
    """Load top-10 from each method - the PUBLISHED shortlists.

    - fixed:   top10_neural_tfs_prioritized.csv (composite + bonuses +
               dual-track 5+5; the real fixed method)
    - centered/uniform: dirichlet_*_top10.csv (composite + bonuses +
               dual-track 5+5)
    All three are the identical quantity class (composite dual-track
    shortlists) - no raw integrated_score fallbacks.
    """
    methods = {}

    fixed_path = RESULTS_DIR / "top10_neural_tfs_prioritized.csv"
    if fixed_path.exists():
        df = pd.read_csv(fixed_path)
        gene_col = "gene_id_v6" if "gene_id_v6" in df.columns else "gene_id"
        methods["fixed"] = set(df[gene_col].astype(str))
    else:
        print(f"WARNING: {fixed_path} missing - the fixed-method arm is "
              f"skipped (run scripts/prioritize_neural_tfs.py first)")

    centered_path = RESULTS_DIR / "dirichlet_centered_top10.csv"
    if centered_path.exists():
        df = pd.read_csv(centered_path)
        gene_col = "gene_id_v6" if "gene_id_v6" in df.columns else "gene_id"
        methods["centered"] = set(df[gene_col].astype(str))

    uniform_path = RESULTS_DIR / "dirichlet_uniform_top10.csv"
    if uniform_path.exists():
        df = pd.read_csv(uniform_path)
        gene_col = "gene_id_v6" if "gene_id_v6" in df.columns else "gene_id"
        methods["uniform"] = set(df[gene_col].astype(str))

    return methods


def universe_size() -> int:
    """Candidate universe N: one row per gene in rank.csv (the shared
    universe of all three methods)."""
    rank = pd.read_csv(RUN_DIR / "rank.csv")
    return int(rank["gene_id"].nunique())


def main():
    print("=== Cross-Method Consensus (valid randomization null) ===")

    methods = load_method_top10()
    print(f"Methods loaded: {list(methods.keys())} (each a dual-track 5+5 shortlist)")

    if len(methods) < 2:
        print("Error: fewer than 2 methods found")
        return 1

    N = universe_size()
    n_methods = len(methods)
    p0 = N_TOP / N
    print(f"Universe N = {N}; chance of landing in one method's top-{N_TOP}: "
          f"p0 = {p0:.3e}")

    all_genes = set()
    for s in methods.values():
        all_genes.update(s)
    all_genes = sorted(all_genes)
    n_genes = len(all_genes)

    # ---- Global set-level overlap (primary statistic) -----------------
    # Pairwise hypergeometric overlap tests between the methods' top-10
    # sets: P(X >= k) under random 10-subsets of the N-gene universe.
    pair_names = []
    pair_overlaps = []
    pair_pvals = []
    method_list = sorted(methods)
    for i in range(len(method_list)):
        for j in range(i + 1, len(method_list)):
            a, b = method_list[i], method_list[j]
            k = len(methods[a] & methods[b])
            p = float(stats.hypergeom.sf(k - 1, N, N_TOP, N_TOP))
            pair_names.append(f"{a}~{b}")
            pair_overlaps.append(k)
            pair_pvals.append(p)
            print(f"  overlap {a} vs {b}: {k}/{N_TOP}  hypergeom p = {p:.3e}")

    # ---- Per-gene consensus strength -----------------------------------
    # Binomial(k successes of n_methods trials at p0 = 10/N). Positively
    # correlated memberships (shared matrix/bonus layer by design) make
    # these descriptive; the set-level hypergeometric is the primary stat.
    consensus_data = []
    for gene in all_genes:
        methods_present = [m for m in methods if gene in methods[m]]
        k = len(methods_present)
        binom_p = _binom_p(k, n_methods, p0)
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

    # Deterministic ordering: k desc -> p_binom asc -> gene_id
    df = df.sort_values(
        ["n_methods_present", "p_binom", "gene_id"],
        ascending=[False, True, True],
    ).reset_index(drop=True)
    df.index = df.index + 1
    df.index.name = "rank"

    print(f"\nTotal genes in any method shortlist: {n_genes}")
    print(f"Genes in >= 2 methods: {(df['n_methods_present'] >= 2).sum()}")
    print(f"Genes in all {n_methods} methods: {(df['n_methods_present'] == n_methods).sum()}")
    print(f"Per-gene significance after Bonferroni (p<0.05): {df['significant_bonferroni'].sum()}")
    print(f"Per-gene significance after BH-FDR (p<0.05): {df['significant_fdr'].sum()}")
    print("(set-level hypergeometric overlaps above are the primary "
          "consensus statistic)")

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
        "universe_size": N,
        "null_p0": p0,
        "n_genes_any_method": n_genes,
        "n_consensus_2plus": int((df["n_methods_present"] >= 2).sum()),
        "n_consensus_all": int((df["n_methods_present"] == n_methods).sum()),
        "n_significant_bonferroni": int(df["significant_bonferroni"].sum()),
        "n_significant_fdr": int(df["significant_fdr"].sum()),
        "pairwise_overlap": {
            name: {"overlap": int(k), "hypergeom_p": p}
            for name, k, p in zip(pair_names, pair_overlaps, pair_pvals)
        },
        "caveat": (
            "The three methods share the candidate matrix and bonus layer "
            "by design, so shortlist memberships are positively correlated; "
            "per-gene binomial p-values are descriptive of agreement "
            "strength. The pairwise hypergeometric set-level overlaps are "
            "the primary consensus statistics."
        ),
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
