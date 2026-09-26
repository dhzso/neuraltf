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
   N = 11,696, p0 = 10/11,696 = 8.5e-4) - NOT the previous binomial
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


def _intersection_pmf(N_s: int, sizes: list[int]) -> np.ndarray:
    """Exact PMF of |A1 ∩ ... ∩ Ak| for k independent uniformly-drawn
    subsets of the given sizes from an N_s-element stratum (2026-09-26).

    Conditioned on its size j, the intersection of uniform subsets is
    itself a uniform j-subset, so intersecting it with a further uniform
    s-subset is Hypergeom(N_s, j, s)."""
    sizes = [int(x) for x in sizes]
    if N_s <= 0 or any(s < 0 or s > N_s for s in sizes):
        raise ValueError(f"invalid stratum draw: N_s={N_s}, sizes={sizes}")
    if any(s == 0 for s in sizes):
        return np.array([1.0])          # intersection deterministically empty
    pmf = np.zeros(sizes[0] + 1)
    pmf[sizes[0]] = 1.0
    for s in sizes[1:]:
        nxt = np.zeros_like(pmf)        # support cannot grow
        for j, w in enumerate(pmf):
            if w == 0.0:
                continue
            ks = np.arange(0, j + 1)
            nxt[:j + 1] += w * stats.hypergeom.pmf(ks, N_s, j, s)
        pmf = nxt
    return pmf


def stratified_overlap_exact_p(k: int, sizes_a: dict[str, int],
                               sizes_b: dict[str, int],
                               N_strata: dict[str, int]) -> float:
    """EXACT one-sided p for the pairwise shortlist overlap under the
    dual-track within-stratum null: the overlap is the sum of
    independent per-stratum Hypergeom(N_s, n_a_s, n_b_s) counts, whose
    exact distribution is a short convolution across strata (2026-09-26;
    replaces the Poisson mean-matched approximation)."""
    dist = np.array([1.0])
    for s, N_s in N_strata.items():
        dist = np.convolve(
            dist, _intersection_pmf(int(N_s), [sizes_a.get(s, 0),
                                               sizes_b.get(s, 0)]))
    if k <= 0:
        return 1.0
    if k >= len(dist):
        return 0.0
    # Sum the UPPER TAIL directly (1 - cumsum catastrophically cancels
    # at extreme k: the observed 10/10 overlaps have true p ~ 1e-26).
    return float(np.clip(dist[k:].sum(), 0.0, 1.0))


def load_method_top10():
    """Load top-10 from each method - the PUBLISHED shortlists.

    - fixed:   top10_neural_tfs_prioritized.csv (composite + bonuses +
               dual-track 5+5; the real fixed method)
    - centered/uniform: dirichlet_*_top10.csv (composite + bonuses +
               dual-track 5+5)
    All three are the identical quantity class (composite dual-track
    shortlists) - no raw integrated_score fallbacks.

    Returns ``(methods, method_sizes)``: gene sets plus each method's
    ACTUAL per-stratum draw sizes (validated in main; 2026-09-26).
    """
    methods = {}
    method_sizes = {}

    def _load(path: Path, name: str):
        if not path.exists():
            if name == "fixed":
                print(f"WARNING: {path} missing - the fixed-method arm is "
                      f"skipped (run scripts/prioritize_neural_tfs.py first)")
            return
        df = pd.read_csv(path)
        gene_col = "gene_id_v6" if "gene_id_v6" in df.columns else "gene_id"
        genes = set(df[gene_col].astype(str))
        methods[name] = genes
        sizes: dict[str, int] = {}
        if "proof_status" in df.columns:
            tested = df["proof_status"].astype(str).eq("tested")
            sizes["tested"] = int(tested.sum())
            sizes["not_tested"] = int((~tested).sum())
        method_sizes[name] = sizes

    _load(RESULTS_DIR / "top10_neural_tfs_prioritized.csv", "fixed")
    _load(RESULTS_DIR / "dirichlet_centered_top10.csv", "centered")
    _load(RESULTS_DIR / "dirichlet_uniform_top10.csv", "uniform")
    return methods, method_sizes


def universe_size() -> int:
    """Candidate universe N: one row per gene in rank.csv (the shared
    universe of all three methods)."""
    rank = pd.read_csv(RUN_DIR / "rank.csv")
    return int(rank["gene_id"].nunique())


def strata_sizes() -> dict[str, int]:
    """Dual-track stratum sizes from rank.csv (2026-09-13 fix; corrected
    2026-09-26).

    The three shortlists are stratified 5 tested + 5 not-tested, not
    uniform random 10-subsets. p0 = 10/N understated per-gene membership
    probability for tested genes by ~73x (5/68 vs 10/11,696) and the
    set-level hypergeometric overstated significance by ~25 orders of
    magnitude. The stratified null is reported alongside the legacy one.

    2026-09-26 fix: the Track-B selection pool is proof_status ==
    "not_tested" EXACTLY (prioritize.assign_tracks puts only these two
    statuses into tracks). The previous len(rank) - n_tested lumped the
    known_fstf genes into the not_tested stratum, understating the
    per-gene null probability by ~0.5%.
    """
    rank = pd.read_csv(RUN_DIR / "rank.csv").drop_duplicates(
        subset="gene_id", keep="first")
    if "proof_status" not in rank.columns:
        return {}
    status = rank["proof_status"].astype(str)
    return {
        "tested": int((status == "tested").sum()),
        "not_tested": int((status == "not_tested").sum()),
    }


def main():
    print("=== Cross-Method Consensus (valid randomization null) ===")

    methods, method_sizes = load_method_top10()
    print(f"Methods loaded: {list(methods.keys())} (each a dual-track 5+5 shortlist)")

    if len(methods) < 2:
        print("Error: fewer than 2 methods found")
        return 1

    # 2026-09-26 validation: confirm the assumed 10-distinct-gene 5+5
    # design against the ACTUAL shortlists (a silent design change would
    # otherwise produce p-values under a wrong design).
    for m in methods:
        n_distinct = len(methods[m])
        sizes = method_sizes.get(m, {})
        if n_distinct != N_TOP or (sizes and set(sizes.values()) != {N_TOP // 2}):
            print(f"WARNING: method '{m}' shortlist is {n_distinct} genes "
                  f"with strata {sizes} (expected 5+5); the nulls below "
                  f"use the ACTUAL sizes")

    N = universe_size()
    n_methods = len(methods)
    p0 = N_TOP / N
    N_strata = strata_sizes()
    # per-stratum membership probability under the dual-track 5+5 design
    # 2026-09-26 fix: fail LOUDLY on an empty stratum — max(N_s, 1)
    # silently produced p0 = 5.0 (a probability > 1).
    p0_strata: dict[str, float] = {}
    if N_strata:
        for s, n_s in N_strata.items():
            if n_s <= 0:
                raise ValueError(
                    f"stratum '{s}' is empty (N_strata={N_strata}); the "
                    "5+5 dual-track null is undefined — check rank.csv "
                    "proof_status values")
            p0_strata[s] = 5 / n_s
    print(f"Universe N = {N}; chance of landing in one method's top-{N_TOP}: "
          f"p0 = {p0:.3e} (legacy uniform null)")
    if p0_strata:
        print(f"Stratified null (dual-track 5+5 design): per-gene p0 = "
              f"{p0_strata} from strata {N_strata}")

    all_genes = set()
    for s in methods.values():
        all_genes.update(s)
    all_genes = sorted(all_genes)
    n_genes = len(all_genes)

    # ---- Global set-level overlap (primary statistic) -----------------
    # Pairwise hypergeometric overlap tests between the methods' top-10
    # sets. THREE nulls are reported (2026-09-26):
    #   legacy:      random 10-subsets of the N-gene universe (upper bound
    #                on significance; ignores the 5+5 stratification)
    #   stratified exact: convolution of the per-stratum hypergeometric
    #                intersection distributions (headline; E[k] =
    #                5^2/n_tested + 5^2/n_not_tested ~ 0.37)
    #   stratified poisson: legacy mean-matched Poisson tail (retained
    #                for transparency only)
    pair_names = []
    pair_overlaps = []
    pair_pvals = []
    pair_pvals_strat = []
    pair_pvals_strat_pois = []
    method_list = sorted(methods)
    lam_strat = None
    if p0_strata:
        lam_strat = sum((5 ** 2) / N_strata[s] for s in N_strata)
    for i in range(len(method_list)):
        for j in range(i + 1, len(method_list)):
            a, b = method_list[i], method_list[j]
            k = len(methods[a] & methods[b])
            # legacy uniform-subset null at the ACTUAL list sizes
            p = float(stats.hypergeom.sf(
                k - 1, N, len(methods[b]), len(methods[a])))
            p_strat = p_strat_pois = None
            if p0_strata:
                sizes_a = method_sizes.get(a, {})
                sizes_b = method_sizes.get(b, {})
                if set(sizes_a) == set(N_strata) and set(sizes_b) == set(N_strata):
                    p_strat = stratified_overlap_exact_p(
                        k, sizes_a, sizes_b, N_strata)
                    p_strat_pois = float(stats.poisson.sf(k - 1, lam_strat))
            pair_names.append(f"{a}~{b}")
            pair_overlaps.append(k)
            pair_pvals.append(p)
            pair_pvals_strat.append(p_strat)
            pair_pvals_strat_pois.append(p_strat_pois)
            msg = f"  overlap {a} vs {b}: {k}/{N_TOP}  hypergeom p = {p:.3e}"
            if p_strat is not None:
                msg += f"  | stratified exact p = {p_strat:.3e}"
            print(msg)

    # ---- Per-gene consensus strength -----------------------------------
    # Binomial(k successes of n_methods trials). Under the stratified
    # design a tested gene's membership probability is 5/n_tested (not
    # 10/N); the per-gene p is computed with the stratum-appropriate p0
    # when the gene's stratum is known from rank.csv. Positively
    # correlated memberships (shared matrix/bonus layer by design) make
    # these descriptive; the set-level hypergeometric is the primary stat.
    rank = pd.read_csv(RUN_DIR / "rank.csv").drop_duplicates(
        subset="gene_id", keep="first")
    # 2026-09-26 fix: map by the EXACT track status. assign_tracks puts
    # only "tested"/"not_tested" genes into tracks, so a shortlist gene
    # with any other status (e.g. known_fstf) means the inputs are
    # inconsistent — fail loudly instead of silently charging it the
    # not_tested p0.
    status_of = {}
    if "proof_status" in rank.columns:
        status_of = dict(zip(
            rank["gene_id"].astype(str),
            rank["proof_status"].astype(str),
        ))
    consensus_data = []
    for gene in all_genes:
        methods_present = [m for m in methods if gene in methods[m]]
        k = len(methods_present)
        p_gene = p0
        if p0_strata:
            st = status_of.get(gene)
            if st in p0_strata:
                p_gene = p0_strata[st]
            elif st is not None:
                raise ValueError(
                    f"shortlist gene {gene} has non-track proof_status "
                    f"{st!r}; the 5+5 dual-track null cannot assign it a "
                    "stratum — check shortlist/rank.csv consistency")
        binom_p = _binom_p(k, n_methods, p_gene)
        consensus_data.append({
            "gene_id": gene,
            "n_methods_present": k,
            "methods": ",".join(sorted(methods_present)),
            "p_binom": binom_p,
            "p0_used": p_gene,
            "is_consensus": k >= 2,
        })

    df = pd.DataFrame(consensus_data)
    pvals = df["p_binom"].values

    bonf_p = bonferroni_correction(pvals)
    fdr_p = benjamini_hochberg(pvals)

    df["p_bonferroni"] = bonf_p
    df["p_fdr_bh"] = fdr_p
    # 2026-09-06 audit: the per-gene significance FLAGS are exported only
    # for genes in >= 2 methods. A single-method gene is a member of the
    # tested family only because it appeared in some shortlist (the family
    # is selected on the outcome), and its binomial p is descriptive of
    # agreement strength, not evidence of consensus. The prior export
    # flagged single-method genes "significant" (e.g. p=0.0026 -> Bonf
    # 0.033) — overclaiming from a self-acknowledged descriptive
    # statistic. Set-level hypergeometric overlaps (above) remain the
    # primary consensus statistics.
    df["significant_bonferroni"] = (bonf_p < 0.05) & (df["n_methods_present"] >= 2)
    df["significant_fdr"] = (fdr_p < 0.05) & (df["n_methods_present"] >= 2)
    df.loc[df["n_methods_present"] < 2, ["p_bonferroni", "p_fdr_bh"]] = np.nan

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
    print(f"Per-gene significance after Bonferroni (p<0.05, k>=2 only): {df['significant_bonferroni'].sum()}")
    print(f"Per-gene significance after BH-FDR (p<0.05, k>=2 only): {df['significant_fdr'].sum()}")
    print("(set-level hypergeometric overlaps above are the primary "
          "consensus statistic; per-gene binomial p-values are descriptive "
          "of agreement strength — the three methods share the candidate "
          "matrix and bonus layer by design)")

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
            name: {"overlap": int(k), "hypergeom_p": p,
                    "hypergeom_p_note": "legacy uniform-subset null (upper bound)",
                    "stratified_exact_p": (float(pe) if pe is not None else None),
                    "stratified_poisson_p": (float(pp) if pp is not None else None)}
            for name, k, p, pe, pp in zip(pair_names, pair_overlaps,
                                          pair_pvals, pair_pvals_strat,
                                          pair_pvals_strat_pois)
        },
        "stratified_null": {
            "design": "dual-track 5 tested + 5 not-tested per method",
            "N_strata": N_strata,
            "p0_strata": p0_strata,
            "expected_pairwise_overlap": (float(lam_strat)
                                           if lam_strat is not None else None),
            "method": (
                "stratified_exact_p: exact convolution of per-stratum "
                "hypergeometric intersection distributions (2026-09-26); "
                "stratified_poisson_p is the legacy mean-matched "
                "approximation, retained for transparency"),
            "note": (
                "The shortlists are stratified (5+5); the stratified null "
                "samples within stratum. Legacy uniform-subset p-values "
                "(hypergeom at 10/N) overstate significance by ~25 orders "
                "of magnitude and are retained only for transparency."),
        },
        "caveat": (
            "The three methods share rank.csv, apply_bonuses(), and "
            "gate_track_b() by design, so shortlist memberships are "
            "positively correlated — the independence assumption of the "
            "binomial and hypergeometric null models is violated. "
            "Per-gene binomial p-values are descriptive of agreement "
            "strength only. All reported p-values are UPPER BOUNDS on "
            "significance (true p-values are larger). The pairwise "
            "hypergeometric set-level overlaps are the primary consensus "
            "statistics. [2026-09-07 audit fix CRITICAL-3]"
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
