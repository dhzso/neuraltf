#!/usr/bin/env python
"""Overlap significance tests between ranking methods.

Tests whether the top-10 overlap between fixed, centered, and uniform
Dirichlet methods is greater than expected by chance using:
- Hypergeometric test for pairwise overlaps (primary statistic)
- Binomial pooled-rate test at the correct null p0 = 10/N
- Three-way overlap with a proper three-set null (NOT the pairwise null
  the previous version used, which was conservative by ~8 orders of
  magnitude for small overlaps)

2026-09-04 alignment: the "fixed" arm now loads the published fixed-method
shortlist (top10_neural_tfs_prioritized.csv) so all three arms are the
same quantity class (composite dual-track shortlists), matching
cross_method_correction.py.

2026-09-26 audit fixes:
- The stratified nulls now report the EXACT tail via convolution of the
  per-stratum hypergeometric intersection distributions (previously a
  Poisson tail with the right mean but the wrong variance — mildly
  anti-conservative at k=1, conservative for large k). The Poisson
  values are retained under *_poisson keys as legacy transparency.
- The 5+5 dual-track design vector is no longer hard-coded: stratum
  draw sizes are derived from each method's ACTUAL shortlist
  composition, and every shortlist is validated (distinct-gene count,
  stratum split), so a silently changed design can no longer produce
  p-values under a wrong design.
- The not_tested stratum size counts proof_status == "not_tested"
  EXACTLY (the rank.csv universe also carries known_fstf genes, which
  are in NEITHER selection track; lumping them into not_tested
  understated the null by ~0.5%).

Usage:
    python scripts/stats/overlap_significance.py
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


def hypergeometric_test(k, n, K, N):
    """P-value for observing >= k successes in sample of n from population of N with K successes."""
    p = stats.hypergeom.sf(k - 1, N, K, n)
    return p


def _strata_of(df: pd.DataFrame, gene_col: str) -> dict[str, set[str]]:
    """Split a dual-track shortlist into its {track: gene set} strata.

    2026-09-13 fix: the three method shortlists are STRATIFIED 5 tested
    (Track A) + 5 not-tested (Track B), not uniform random 10-subsets of
    the ~11.7k-gene universe. The old null (10 random genes from N)
    understated the expected overlap by ~40x (E[overlap] ~ 0.0086 under
    the old null vs ~0.37 under the correct within-stratum null for
    tested genes: 5*5/n_tested + 5*5/n_not_tested), overstating
    significance by ~25-38 orders of magnitude. The corrected null
    samples within stratum, independently per method.
    """
    strata = {}
    if "proof_status" in df.columns:
        tested = df["proof_status"].astype(str).eq("tested")
        strata["tested"] = set(df.loc[tested, gene_col].astype(str))
        strata["not_tested"] = set(df.loc[~tested, gene_col].astype(str))
    else:
        strata["all"] = set(df[gene_col].astype(str))
    return strata


def three_way_null_p(k3, n, N, n_methods=3):
    """Exact null p for k3 genes shared by all n_methods independent
    random n-subsets of an N-gene universe.

    P(a given gene is in one random n-subset) = n/N; independence across
    methods gives P(in all n_methods) = (n/N)^n_methods, and the count
    of 3-way-shared genes is Binomial(N, (n/N)^n_methods) under the
    independent-random-subsets null. One-sided P(X >= k3).

    2026-09-13: superseded for the dual-track shortlists by
    ``three_way_stratified_null_p`` (kept for the uniform-subset case).
    """
    p_each = (n / N) ** n_methods
    if hasattr(stats, "binomtest"):
        return float(stats.binomtest(k3, N, p_each, alternative="greater").pvalue)
    return float(stats.binom_test(k3, N, p_each, alternative="greater"))


def expected_stratified_overlap(n_a: dict[str, int], n_b: dict[str, int],
                                N_strata: dict[str, int]) -> float:
    """Expected pairwise overlap under the dual-track within-stratum null.

    Each method draws its n_s genes uniformly from stratum s (sizes
    N_s); the expected overlap of the two independent draws is
    sum_s n_a[s] * n_b[s] / N_s.
    """
    return sum(
        (n_a.get(s, 0) * n_b.get(s, 0)) / max(N_strata[s], 1)
        for s in N_strata
    )


def _intersection_pmf(N_s: int, sizes: list[int]) -> np.ndarray:
    """Exact PMF of |A1 ∩ ... ∩ Ak| for k independent uniformly-drawn
    subsets of the given sizes from an N_s-element stratum.

    Conditioned on its size j, the intersection of uniform subsets is
    itself a uniform j-subset, so intersecting it with a further uniform
    s-subset is Hypergeom(N_s, j, s). Returns an array indexed by
    intersection size (0..max).
    """
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


def _stratified_total_pmf(N_strata: dict[str, int],
                          sizes_by_method: list[dict[str, int]]) -> np.ndarray:
    """PMF of the total overlap count (sum across strata) under the
    within-stratum null: convolution of the per-stratum intersection
    PMFs across independent strata."""
    dist = np.array([1.0])
    for s, N_s in N_strata.items():
        sizes = [m.get(s, 0) for m in sizes_by_method]
        dist = np.convolve(dist, _intersection_pmf(int(N_s), sizes))
    return dist


def _tail_p(dist: np.ndarray, k: int) -> float:
    """One-sided P(X >= k) for a discrete PMF array.

    Sums the UPPER TAIL directly: the complement form
    ``1 - dist[:k].sum()`` catastrophically cancels for extreme k (the
    observed 10/10 method overlaps have true p ~ 1e-26, which the
    complement rounded to the ~1e-16 floating-point floor).
    """
    if k <= 0:
        return 1.0
    if k >= len(dist):
        return 0.0
    return float(np.clip(dist[k:].sum(), 0.0, 1.0))


def stratified_overlap_p(k_observed: int, sizes_a: dict[str, int],
                         sizes_b: dict[str, int],
                         N_strata: dict[str, int]) -> float:
    """EXACT one-sided p for the pairwise overlap under the dual-track
    within-stratum null (2026-09-26): the overlap count is the sum of
    independent per-stratum Hypergeom(N_s, n_a_s, n_b_s) counts, whose
    exact distribution is a short convolution. Replaces the Poisson
    approximation, which had the correct mean but the wrong variance.
    """
    return _tail_p(_stratified_total_pmf(N_strata, [sizes_a, sizes_b]),
                   k_observed)


def stratified_overlap_p_poisson(k_observed: int, sizes_a: dict[str, int],
                                sizes_b: dict[str, int],
                                N_strata: dict[str, int]) -> float:
    """LEGACY Poisson approximation of the stratified pairwise-overlap
    tail (retained for transparency; superseded by the exact
    convolution in ``stratified_overlap_p``)."""
    lam = expected_stratified_overlap(sizes_a, sizes_b, N_strata)
    if lam <= 0:
        return 1.0
    return float(stats.poisson.sf(k_observed - 1, lam))


def three_way_stratified_null_p(k3: int, sizes_a: dict[str, int],
                                sizes_b: dict[str, int],
                                sizes_c: dict[str, int],
                                N_strata: dict[str, int]) -> float:
    """EXACT one-sided p for the 3-way overlap under the dual-track
    within-stratum null (2026-09-26): per stratum the 3-way-shared count
    is a hypergeometric chain (|A∩B| ~ Hypergeom(N_s, a, b), then
    |A∩B∩C| | |A∩B|=j ~ Hypergeom(N_s, j, c)); the total across strata
    is the convolution. Replaces the Poisson approximation."""
    return _tail_p(
        _stratified_total_pmf(N_strata, [sizes_a, sizes_b, sizes_c]), k3)


def three_way_stratified_null_p_poisson(k3: int, sizes_a: dict[str, int],
                                        sizes_b: dict[str, int],
                                        sizes_c: dict[str, int],
                                        N_strata: dict[str, int]) -> float:
    """LEGACY Poisson approximation of the 3-way stratified tail
    (retained for transparency; superseded by the exact convolution in
    ``three_way_stratified_null_p``)."""
    lam = sum(
        (sizes_a.get(s, 0) * sizes_b.get(s, 0) * sizes_c.get(s, 0))
        / max(N_strata[s], 1) ** 2
        for s in N_strata
    )
    if lam <= 0:
        return 1.0
    return float(stats.poisson.sf(k3 - 1, lam))


def _binom_p(k, n, p):
    """Compute binomial test p-value compatible with SciPy <1.12 and >=1.12."""
    if hasattr(stats, "binomtest"):
        return float(stats.binomtest(k, n, p, alternative="greater").pvalue)
    return float(stats.binom_test(k, n, p, alternative="greater"))


def main():
    print("=== Overlap Significance Tests ===")

    centered_path = RESULTS_DIR / "dirichlet_centered_top10.csv"
    uniform_path = RESULTS_DIR / "dirichlet_uniform_top10.csv"
    fixed_path = RESULTS_DIR / "top10_neural_tfs_prioritized.csv"

    methods = {}
    method_strata = {}

    def _gene_col(df):
        return "gene_id_v6" if "gene_id_v6" in df.columns else (
            "gene_id" if "gene_id" in df.columns else "gene_id_v6")

    if centered_path.exists():
        df_c = pd.read_csv(centered_path)
        methods["centered"] = set(df_c[_gene_col(df_c)].astype(str))
        method_strata["centered"] = _strata_of(df_c, _gene_col(df_c))
    else:
        print(f"Warning: {centered_path} not found")
        methods["centered"] = set()
        method_strata["centered"] = {}

    if uniform_path.exists():
        df_u = pd.read_csv(uniform_path)
        methods["uniform"] = set(df_u[_gene_col(df_u)].astype(str))
        method_strata["uniform"] = _strata_of(df_u, _gene_col(df_u))
    else:
        print(f"Warning: {uniform_path} not found")
        methods["uniform"] = set()
        method_strata["uniform"] = {}

    if fixed_path.exists():
        # The PUBLISHED fixed-method shortlist (composite dual-track),
        # not a raw integrated_score top-10 from rank.csv.
        df_f = pd.read_csv(fixed_path)
        methods["fixed"] = set(df_f[_gene_col(df_f)].astype(str))
        method_strata["fixed"] = _strata_of(df_f, _gene_col(df_f))
    else:
        print(f"Warning: {fixed_path} not found")
        methods["fixed"] = set()
        method_strata["fixed"] = {}

    # Population size: the SHARED candidate universe - one row per gene in
    # rank.csv (never the row-exploded Dirichlet CSVs).
    rank_csv = RUN_DIR / "rank.csv"
    rank_df = pd.read_csv(rank_csv) if rank_csv.exists() else pd.DataFrame(
        {"gene_id": [], "proof_status": []})
    # 2026-09-26: aligned with power_analysis.py's emergency fallback.
    N = rank_df["gene_id"].nunique() if len(rank_df) else 11696
    n = 10

    # 2026-09-13 stratified null: stratum sizes from rank.csv itself (the
    # label the shortlists stratify on), not guessed constants.
    # 2026-09-26 fix: the Track-B pool is proof_status == "not_tested"
    # EXACTLY — the rank.csv universe also carries known_fstf genes,
    # which sit in NEITHER selection track (prioritize.assign_tracks
    # selects only tested / not_tested). Lumping them into not_tested
    # understated the null draw probability by ~0.5%.
    if "proof_status" in rank_df.columns and len(rank_df):
        status = rank_df["proof_status"].astype(str)
        n_tested = int((status == "tested").sum())
        n_not_tested = int((status == "not_tested").sum())
        N_strata = {"tested": n_tested, "not_tested": n_not_tested}
        for s, n_s in N_strata.items():
            if n_s <= 0:
                raise ValueError(
                    f"stratum '{s}' is empty in rank.csv (N_strata={N_strata}); "
                    "the dual-track within-stratum null is undefined — "
                    "check rank.csv proof_status values")
    else:
        N_strata = {"all": N}

    # 2026-09-26 fix: derive the per-method draw sizes from the ACTUAL
    # shortlist composition instead of assuming the hard-coded 5+5
    # design, and validate each shortlist. A silently changed design
    # (a gate change, a corrupted CSV) previously produced p-values
    # under a wrong design with no diagnostic.
    method_sizes: dict[str, dict[str, int]] = {}
    method_names = list(methods.keys())
    print(f"\nMethods found: {method_names}")
    for m in method_names:
        sizes = {s: len(g) for s, g in method_strata[m].items()} \
            if method_strata.get(m) else {}
        method_sizes[m] = sizes
        n_distinct = len(methods[m])
        design_note = ", ".join(f"{s}={sizes[s]}" for s in sorted(sizes)) \
            or "unstratified"
        if n_distinct != n:
            print(f"  WARNING: method '{m}' shortlist has {n_distinct} "
                  f"distinct genes (expected {n}); the nulls below use "
                  f"the ACTUAL sizes")
        print(f"  {m}: {n_distinct} genes, strata {design_note}")
    # Expected pairwise overlap for a representative pair (all pairs
    # share the design when the shortlists conform; uses the first two
    # loaded methods otherwise).
    if "tested" in N_strata:
        ref = [method_sizes[m] for m in method_names[:2]] or []
        exp_pair = (expected_stratified_overlap(ref[0], ref[-1], N_strata)
                    if ref else float("nan"))
        print(f"  Strata universe: {N_strata} "
              f"(expected pairwise overlap under null: {exp_pair:.3f})")

    results = {
        "pairwise": {}, "three_way": {}, "binomial": {}, "overlaps": {},
        "stratified_null": {},
    }
    if "tested" in N_strata:
        ref = [method_sizes[m] for m in method_names[:2]] or []
        results["stratified_null"] = {
            "design": "dual-track 5 tested + 5 not-tested per method",
            "N_strata": N_strata,
            "n_strata_per_method": method_sizes,
            "expected_pairwise_overlap": (
                expected_stratified_overlap(ref[0], ref[-1], N_strata)
                if ref else None),
            "method": (
                "exact: convolution of per-stratum hypergeometric "
                "intersection distributions (2026-09-26); poisson keys "
                "are the legacy mean-matched approximation"),
            "note": (
                "The shortlists are stratified (5+5), so the correct null "
                "samples within stratum. The legacy uniform-subset null "
                "(10 random of N) overstated significance by ~25-38 orders "
                "of magnitude; both are reported for transparency."),
        }

    print("\n--- Pairwise Overlaps ---")
    for i in range(len(method_names)):
        for j in range(len(method_names)):
            if i == j:
                continue
            m1, m2 = method_names[i], method_names[j]
            overlap = methods[m1] & methods[m2]
            union = methods[m1] | methods[m2]
            k = len(overlap)
            jaccard = float(k / len(union)) if union else 0.0
            # 2026-09-26: legacy uniform-subset null parameterized by the
            # ACTUAL list sizes (equals the old fixed-n form at n=10).
            hg_p = hypergeometric_test(
                k, len(methods[m1]), len(methods[m2]), N)
            # 2026-09-26: stratification-corrected p — EXACT convolution
            # of the per-stratum hypergeometric intersection null
            # (headline) plus the legacy Poisson tail, both parameterized
            # by each method's ACTUAL stratum draw sizes.
            sizes_a = method_sizes.get(m1, {})
            sizes_b = method_sizes.get(m2, {})
            strat_keys_ok = (
                "tested" in N_strata
                and set(sizes_a) == set(N_strata)
                and set(sizes_b) == set(N_strata)
            )
            strat_p = stratified_overlap_p(k, sizes_a, sizes_b, N_strata) \
                if strat_keys_ok else None
            strat_p_pois = stratified_overlap_p_poisson(
                k, sizes_a, sizes_b, N_strata) if strat_keys_ok else None

            key = f"{m1}_vs_{m2}"
            results["pairwise"][key] = {
                "overlap_count": k,
                "jaccard": jaccard,
                "overlap_genes": sorted(list(overlap)),
                "hypergeometric_p": float(hg_p),
                "hypergeometric_p_note": "legacy uniform-subset null (upper bound)",
                "stratified_exact_p": strat_p,
                "stratified_poisson_p": strat_p_pois,
                "N_population": N,
            }
            if i < j:
                msg = (f"  {m1} vs {m2}: {k}/10 overlap (Jaccard={jaccard:.2f}), "
                       f"legacy hypergeom p={hg_p:.4e}")
                if strat_p is not None:
                    msg += f", stratified exact p={strat_p:.4e}"
                print(msg)

    for k, v in results["pairwise"].items():
        results["overlaps"][k] = {
            "count": v["overlap_count"],
            "jaccard": v["jaccard"],
            "p_value": v["hypergeometric_p"],
        }

    print("\n--- Three-way Overlap ---")
    if len(methods) == 3:
        three_way = methods.get("centered", set()) & methods.get("uniform", set()) & methods.get("fixed", set())
        k3 = len(three_way)
        # 2026-09-26: dual-track corrected 3-set null — EXACT convolution
        # of the per-stratum hypergeometric chains (headline) plus the
        # legacy Poisson tail; both parameterized by each method's
        # ACTUAL stratum draw sizes. The legacy uniform-subset
        # Binomial(N, (n/N)^3) is kept for transparency and labelled as
        # an upper bound on significance.
        sizes_3way = [method_sizes.get(m, {}) for m in ("centered", "uniform", "fixed")]
        ok_3way = (
            "tested" in N_strata
            and all(set(x) == set(N_strata) for x in sizes_3way)
        )
        p3_strat = three_way_stratified_null_p(
            k3, sizes_3way[0], sizes_3way[1], sizes_3way[2], N_strata
        ) if ok_3way else None
        p3_pois = three_way_stratified_null_p_poisson(
            k3, sizes_3way[0], sizes_3way[1], sizes_3way[2], N_strata
        ) if ok_3way else None
        p3_legacy = three_way_null_p(k3, n, N)
        results["three_way"] = {
            "overlap_count": k3,
            "overlap_genes": sorted(list(three_way)),
            "legacy_null_p_per_gene": float((n / N) ** 3),
            "binomial_three_set_p": float(p3_legacy),
            "stratified_exact_p": p3_strat,
            "stratified_poisson_p": p3_pois,
            "expected_stratified_overlap": (
                float(sum(
                    sizes_3way[0].get(s, 0) * sizes_3way[1].get(s, 0)
                    * sizes_3way[2].get(s, 0) / max(N_strata[s], 1) ** 2
                    for s in N_strata)) if ok_3way else None),
            "N_population": N,
        }
        results["overlaps"]["three_way"] = {
            "count": k3,
            "p_value": float(p3_strat if p3_strat is not None else p3_legacy),
        }
        msg = f"  Three-way overlap: {k3}/10, legacy 3-set p={p3_legacy:.4e}"
        if p3_strat is not None:
            msg += f", stratified exact p={p3_strat:.4e}"
        print(msg)

    print("\n--- Binomial Test (legacy uniform-subset null) ---")
    unique_pairs = [(method_names[i], method_names[j]) for i in range(len(method_names)) for j in range(i+1, len(method_names))]
    total_possible_pairs = len(unique_pairs)
    total_overlap_count = sum(
        results["pairwise"][f"{m1}_vs_{m2}"]["overlap_count"]
        for m1, m2 in unique_pairs
        if f"{m1}_vs_{m2}" in results["pairwise"]
    )
    # 2026-09-26: slot budget and null rate from the ACTUAL list sizes
    # (each pair can share at most min(|a|, |b|) genes; the per-slot
    # null rate is the expected overlap fraction under independence).
    max_possible = sum(
        min(len(methods[m1]), len(methods[m2]))
        for m1, m2 in unique_pairs
        if m1 in methods and m2 in methods
    )
    expected_total = sum(
        len(methods[m1]) * len(methods[m2]) / N
        for m1, m2 in unique_pairs
        if m1 in methods and m2 in methods
    )
    if max_possible > 0:
        null_rate = expected_total / max_possible
        binom_p = _binom_p(total_overlap_count, max_possible, null_rate)
        results["binomial"] = {
            "total_overlaps": total_overlap_count,
            "max_possible": max_possible,
            "observed_rate": total_overlap_count / max_possible,
            "expected_rate_under_null": float(null_rate),
            "binomial_p_greater": float(binom_p),
        }
        print(f"  Total overlaps: {total_overlap_count}/{max_possible}, binom p={binom_p:.4e}")
    else:
        results["binomial"] = {"error": "no pairs"}

    # 2026-09-07 audit fix (CRITICAL-3): the three methods share rank.csv,
    # apply_bonuses(), and gate_track_b() — they are NOT independent.
    # The hypergeometric/binomial null understates expected overlap,
    # inflating significance. Flag this in the output for consumers.
    results["CAVEAT_independence_assumption"] = (
        "These p-values assume method independence, which is violated: "
        "all 3 methods share rank.csv, apply_bonuses(), and gate_track_b(). "
        "The null model understates expected overlap. Interpret as an UPPER "
        "BOUND on significance (true p-values are larger)."
    )

    out_path = RESULTS_DIR / "overlap_significance.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
