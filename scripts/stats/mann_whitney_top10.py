#!/usr/bin/env python
"""Mann-Whitney U test for the published Top-10 shortlist vs remaining candidates.

2026-09-04 audit note (selection tautology, documented): a top-10 group
defined as the 10 highest values of the SAME score being tested makes
the U-test tautological - under any score distribution the top-10
exceed the rest, U = n1*n2 and delta = 1 by construction.

2026-09-19 reference-set fix: the "top-10" group was previously the 10
highest INTEGRATED scores (only 5/10 shared with the published
dual-track shortlist). It is now the PUBLISHED shortlist
(results/top10_neural_tfs_prioritized.csv). The circular/integrated arm
remains diagnostic (composite = base + bonuses); the honest/strict arms
contrast the shortlist against a label-free score the shortlist was NOT
selected by, so their separation is genuine evidence.

The script therefore reports three comparisons:
1. circular: shortlist vs rest on integrated_score (near-tautological;
   kept for continuity with figure 23's dual circular/honest framing).
2. honest: shortlist vs rest on a RECOMPUTED score that EXCLUDES the
   streams directly derived from the ground-truth label (rnai,
   neural_enriched, neural_specificity, perez_lineage) - the
   label-leakage-free contrast (same recomputation as
   precision_recall.py).
3. strict: additionally excludes reproducibility (lower bound).

Usage:
    python scripts/stats/mann_whitney_top10.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from bioforge.evidence.scoring import DEFAULT_WEIGHTS, STREAM_ORDER  # noqa: E402

RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
FIG_DIR = REPO / "projects" / "NeuralTF" / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

STREAMS = [getattr(s, "value", s) for s in STREAM_ORDER]
W_BY_NAME = {getattr(k, "value", k): float(v)
             for k, v in DEFAULT_WEIGHTS.items()}
# Label-derived streams excluded from the honest score. 2026-09-06 audit:
# perez_lineage alone has AUC 0.9903 on the RNAi-validated label (it is a
# hand-curated neural-TF-family membership — the same species of leakage
# as rnai), and reproducibility alone has AUC 0.9078 (atlas_membership
# embeds King-neural-G0 hits). The honest score now excludes ALL of
# these; a second, stricter variant also drops reproducibility.
HONEST_EXCLUDE = ["rnai", "neural_enriched", "neural_specificity",
                  "perez_lineage"]
HONEST_STRICT_EXCLUDE = HONEST_EXCLUDE + ["reproducibility"]


def rank_biserial_correlation(u_stat, n1, n2):
    """Rank-biserial correlation from Mann-Whitney U.

    U counts pairs where a group-1 value EXCEEDS a group-2 value (plus
    half the ties). The conventional effect size favoring group 1 is
    r = 2U/(n1*n2) - 1, giving +1 for perfect separation (group 1 all
    higher) and -1 for the reverse.
    """
    return (2.0 * u_stat) / (n1 * n2) - 1.0


def honest_score(df: pd.DataFrame, exclude: list[str] | None = None) -> pd.Series:
    """Renormalized weighted score over label-independent streams only
    (identical to EvidenceScorer's per-record renormalization).

    ``exclude`` defaults to HONEST_EXCLUDE; pass HONEST_STRICT_EXCLUDE
    for the reproducibility-free variant.
    """
    excl = HONEST_EXCLUDE if exclude is None else exclude
    keep = [s for s in STREAMS if s in df.columns and s not in excl]
    W = np.array([W_BY_NAME[s] for s in keep])
    M = df[keep].to_numpy(dtype=float)
    mask = ~np.isnan(M)
    num = np.where(mask, M, 0.0) @ W
    den = (mask * W).sum(axis=1)
    return pd.Series(np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0),
                     index=df.index)


def compare(top_scores, rest_scores):
    n1, n2 = len(top_scores), len(rest_scores)
    u_stat, p_one = stats.mannwhitneyu(top_scores, rest_scores, alternative="greater")
    _, p_two = stats.mannwhitneyu(top_scores, rest_scores, alternative="two-sided")
    r_rb = rank_biserial_correlation(u_stat, n1, n2)
    cliff = 0.0
    for t in top_scores:
        cliff += np.sum(t > rest_scores) - np.sum(t < rest_scores)
    cliff /= (n1 * n2)
    return {
        "top10_n": int(n1), "rest_n": int(n2),
        "top10_mean": float(np.mean(top_scores)),
        "rest_mean": float(np.mean(rest_scores)),
        "top10_median": float(np.median(top_scores)),
        "rest_median": float(np.median(rest_scores)),
        "u_statistic": float(u_stat),
        "p_value_one_sided": float(p_one),
        "p_value_two_sided": float(p_two),
        "rank_biserial_correlation": float(r_rb),
        "cliffs_delta": float(cliff),
    }


def main():
    print("=== Mann-Whitney U: Top-10 vs Rest (circular + honest) ===")

    p = RUN_DIR / "rank.csv"
    if not p.exists():
        print("Error: rank.csv not found")
        return 1
    df = pd.read_csv(p).drop_duplicates(subset="gene_id", keep="first")

    score_col = "integrated_score"
    gene_col = "gene_id"
    df = df.dropna(subset=[score_col])

    # 2026-09-19: the top-10 group is the PUBLISHED dual-track shortlist
    # (not the integrated-score top-10, which shares only 5/10 members).
    shortlist_path = RESULTS_DIR / "top10_neural_tfs_prioritized.csv"
    if not shortlist_path.exists():
        print("Error: run scripts/prioritize_neural_tfs.py first "
              "(top10_neural_tfs_prioritized.csv missing)")
        return 1
    published = pd.read_csv(shortlist_path)
    id_col = "gene_id_v6" if "gene_id_v6" in published.columns else "gene_id"
    top10_genes = set(published[id_col].astype(str))

    # 1) circular comparison (tautological by construction)
    top10_scores = df.loc[df[gene_col].isin(top10_genes), score_col].values
    rest_scores = df.loc[~df[gene_col].isin(top10_genes), score_col].values
    circular = compare(top10_scores, rest_scores)
    print(f"\n[circular] Top-10: n={circular['top10_n']}, "
          f"mean={circular['top10_mean']:.4f}  |  Rest: n={circular['rest_n']}, "
          f"mean={circular['rest_mean']:.4f}")
    print(f"[circular] U={circular['u_statistic']:.0f}  p={circular['p_value_one_sided']:.3e}  "
          f"r_rb={circular['rank_biserial_correlation']:.3f}  delta={circular['cliffs_delta']:.3f}"
          f"  (near-tautological: shortlist composite is dominated by the tested base score)")

    # 2) honest comparison: same shortlist, label-independent score.
    # Two honesty levels (2026-09-06 audit):
    #   honest  — excludes rnai/neural_*/perez_lineage
    #   strict  — additionally excludes reproducibility (King-membership
    #             leakage); the lower bound on true discrimination.
    hs = honest_score(df)
    top10_h = hs[df[gene_col].isin(top10_genes)].values
    rest_h = hs[~df[gene_col].isin(top10_genes)].values
    honest = compare(top10_h, rest_h)
    print(f"\n[honest]   Top-10: mean={honest['top10_mean']:.4f}  |  "
          f"Rest: mean={honest['rest_mean']:.4f}")
    print(f"[honest]   U={honest['u_statistic']:.0f}  p={honest['p_value_one_sided']:.3e}  "
          f"r_rb={honest['rank_biserial_correlation']:.3f}  delta={honest['cliffs_delta']:.3f}"
          f"  (score excludes rnai/neural_*/perez_lineage)")

    ss = honest_score(df, exclude=HONEST_STRICT_EXCLUDE)
    top10_s = ss[df[gene_col].isin(top10_genes)].values
    rest_s = ss[~df[gene_col].isin(top10_genes)].values
    strict = compare(top10_s, rest_s)
    print(f"[strict]   U={strict['u_statistic']:.0f}  p={strict['p_value_one_sided']:.3e}  "
          f"r_rb={strict['rank_biserial_correlation']:.3f}  delta={strict['cliffs_delta']:.3f}"
          f"  (additionally excludes reproducibility)")

    results = {
        "test": "Mann-Whitney U (one-sided, greater)",
        "caveat": (
            "The 'circular' comparison contrasts the PUBLISHED dual-track "
            "shortlist against the integrated base score that dominates "
            "its own composite — DIAGNOSTIC ONLY, near-tautological, "
            "carries zero evidential value for the selection itself. The "
            "'honest' comparison recomputes the score excluding "
            "rnai/neural_enriched/neural_specificity/perez_lineage "
            "(label-derived streams; perez_lineage alone carries AUC 0.99 "
            "on the label) — the shortlist was NOT selected on this "
            "score, so its separation is genuine evidence. The 'strict' "
            "variant additionally excludes reproducibility (King-neural-G0 "
            "membership leakage) and is the lower bound on true "
            "discrimination. [2026-09-19 reference-set fix]"
        ),
        "circular": circular,
        "honest": honest,
        "honest_strict": strict,
        "honest_exclude": HONEST_EXCLUDE,
        "honest_strict_exclude": HONEST_STRICT_EXCLUDE,
        "top10_genes": sorted(list(top10_genes)),
        "label_note": (
            "2026-09-11 ground-truth correction: the historical label "
            "('tested') is the King mmc5 RNAi SCREENING list — mmc5 is "
            "titled 'All Transcription Factors Inhibited' and the "
            "distributed copy lost the red/green phenotype font encoding, "
            "so membership means RNAi was performed, not that a phenotype "
            "was observed."
        ),
    }

    print(f"\nTop-10 genes (published dual-track shortlist):")
    for _, row in published.iterrows():
        gid = str(row[id_col])
        nm = row.get("gene_name", gid)
        nm = nm if isinstance(nm, str) and str(nm) != "nan" else gid
        s = df.loc[df[gene_col] == gid, score_col]
        s_val = float(s.iloc[0]) if len(s) else float("nan")
        print(f"  {str(nm)[:28]:>28}  integrated={s_val:.4f}  "
              f"composite={row['composite_score']:.4f}")

    out_path = RESULTS_DIR / "mann_whitney_top10.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

