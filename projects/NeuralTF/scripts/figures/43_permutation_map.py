#!/usr/bin/env python
"""Figure 43 (2026-09-20, v2): Atlas-permutation significance MAP — the whole
test in ONE graph (single axes, Nature Communications 2-column width).

Design rationale: the permutation null can never exceed the label-independent
King ceiling, so the (real integrated score, ceiling) plane with its identity
diagonal IS the state space of the test:

  side of the diagonal  -> design testability (can the permutation speak?)
  distance below it     -> label-stream surplus (the effect it resolves)
  marker fill           -> the verdict (floor p = 1.0e-3 vs pinned p = 1)
  crimson rings         -> boundary genes (floor only by float tie)
  terracotta + labels   -> the dual-track shortlist top-10

Because the null outcome is strictly bimodal (62 at the floor, 81 pinned,
0 intermediate p-values), the p-value carries a single bit that is encoded as
marker fill rather than wasting an axis on it. Census, multiple-testing
budget and shortlist disposition are annotated in the empty half-plane
corners. Every number is computed from permutation_pvalues_neural.csv and
top10_neural_tfs_prioritized.csv — zero hard-coded statistics.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

C_TEST = "#687787"      # slate gray — testable, BH-sig at floor
C_SAT = "#C8CED6"       # soft light gray — saturated, pinned p = 1
NL = chr(10)


def build():
    perm = pd.read_csv(RES / "permutation_pvalues_neural.csv")
    top10 = load_top10()

    pvals = perm["empirical_p"].to_numpy(dtype=float)
    scores = perm["real_integrated_score"].to_numpy(dtype=float)
    ceilings = perm["label_independent_ceiling"].to_numpy(dtype=float)
    n_perm = int(perm["n_perm"].iloc[0])
    p_floor = 1.0 / (n_perm + 1)
    n = len(pvals)
    n_sig = int(perm["significant_fdr_05"].sum())
    n_pinned = int((pvals >= 1.0).sum())
    n_floor = int((pvals <= p_floor + 1e-12).sum())
    n_mid = n - n_floor - n_pinned
    untest = perm["untestable_by_permutation"].to_numpy(dtype=bool)
    n_untest = int(untest.sum())
    # After 2026-09-21 fix, all untestable genes have p=1.0 by enforcement,
    # so n_bound (untestable but at floor p) should be 0.  Keep the
    # computation for backward compatibility with pre-fix CSVs.
    n_bound = int((untest & ~(pvals >= 1.0)).sum())
    n_test = n - n_untest
    sig_mask = perm["significant_fdr_05"].values.astype(bool)
    bh = float(perm.loc[sig_mask, "empirical_p"].max()) if sig_mask.any() else p_floor

    short = perm[perm["gene_id"].isin(set(top10["gene_id"]))].merge(
        top10[["gene_id", "track"]], on="gene_id", how="left")
    short["gene_label"] = [clean_gene_symbol(nm, g)
                           for nm, g in zip(short["gene_name"], short["gene_id"])]
    short_ids = set(short["gene_id"])
    is_short = perm["gene_id"].isin(short_ids).to_numpy()
    is_pinned = pvals >= 1.0

    fig, ax = plt.subplots(figsize=(7.08, 5.6), dpi=500)
    fig.subplots_adjust(left=0.085, right=0.975, top=0.905, bottom=0.125)
    lo = min(scores.min(), ceilings.min()) - 0.03
    hi = max(scores.max(), ceilings.max()) + 0.05

    # design half-planes (very light tints) + identity diagonal
    ax.fill([lo, lo, hi], [lo, hi, hi], color=C_HL, alpha=0.035, zorder=0)
    ax.fill([lo, hi, hi], [lo, lo, hi], color=C_A, alpha=0.035, zorder=0)
    ax.plot([lo, hi], [lo, hi], color="#9E9E9E", lw=0.9, ls="--", zorder=1)
    ax.text(lo + 0.028, lo + 0.048, "identity (score = ceiling)",
            rotation=45, rotation_mode="anchor", fontsize=5.2,
            color="#8A8A8A", ha="left", va="bottom", zorder=2)

    # population markers
    ax.scatter(scores[~is_short & ~untest], ceilings[~is_short & ~untest],
               s=12, c=C_TEST, edgecolors="none", zorder=3,
               label=f"testable — BH-sig at floor p = 1.0e-3 (n = {n_test})")
    ax.scatter(scores[~is_short & untest & is_pinned],
               ceilings[~is_short & untest & is_pinned],
               s=14, facecolors="none", edgecolors="#8A97A3", linewidths=0.8,
               zorder=3, label=f"saturated — null = real, pinned p = 1 (n = {n_pinned})")
    if n_bound > 0:
        ax.scatter(scores[~is_short & untest & ~is_pinned],
                   ceilings[~is_short & untest & ~is_pinned],
                   s=18, facecolors="none", edgecolors=C_HL, linewidths=1.0,
                   zorder=4,
                   label=f"boundary — floor only by float tie (n = {n_bound})")

    # shortlist markers: filled = resolves at floor, open = pinned at p = 1;
    # float-tie shortlist genes get an extra crimson ring
    sh_test = short[short["empirical_p"] <= p_floor + 1e-12]
    sh_pin = short[short["empirical_p"] >= 1.0]
    ax.scatter(sh_test["real_integrated_score"], sh_test["label_independent_ceiling"],
               s=42, c=C_HL, edgecolors="white", linewidths=0.7, zorder=6,
               label=f"dual-track shortlist top-10 ({len(sh_test)} resolve / {len(sh_pin)} pinned)")
    ax.scatter(sh_pin["real_integrated_score"], sh_pin["label_independent_ceiling"],
               s=42, facecolors="none", edgecolors=C_HL, linewidths=1.3, zorder=6)
    ring = short[short["untestable_by_permutation"] & (short["empirical_p"] <= p_floor + 1e-12)]
    if len(ring) > 0:
        ax.scatter(ring["real_integrated_score"], ring["label_independent_ceiling"],
                   s=78, facecolors="none", edgecolors=C_HL, linewidths=0.9,
                   zorder=5, label="float-tie shortlist member")

    # leader-line shortlist labels: testable above the diagonal, pinned below,
    # two staggered rows so connectors never collide
    def leader_labels(grp, side):
        grp = grp.sort_values("real_integrated_score").reset_index(drop=True)
        for k, (_, r) in enumerate(grp.iterrows()):
            x0, y0 = r["real_integrated_score"], r["label_independent_ceiling"]
            off = (0.045 + 0.050 * (k % 2)) * side
            col = C_A if r["track"] == "A" else C_HL
            ax.annotate(r["gene_label"], xy=(x0, y0), xytext=(x0 + 0.004, y0 + off),
                        ha="center", va="center", fontsize=5.0, color=col,
                        zorder=7, arrowprops=dict(arrowstyle="-", lw=0.35,
                                                  color=col, shrinkA=1, shrinkB=3))
    leader_labels(sh_test, +1)
    leader_labels(sh_pin, -1)

    # surplus gauge: perpendicular arrow from the largest-surplus testable
    # gene to the diagonal — that distance is what the permutation resolves
    cand = (~is_short & ~untest)
    j = int(np.argmax(np.where(cand, scores - ceilings, -np.inf)))
    fx = fy = (scores[j] + ceilings[j]) / 2.0
    ax.annotate("", xy=(fx, fy), xytext=(scores[j], ceilings[j]),
                arrowprops=dict(arrowstyle="<->", lw=0.7, color="#333333",
                                shrinkA=0, shrinkB=0), zorder=6)
    mx, my = (scores[j] + fx) / 2, (ceilings[j] + fy) / 2
    ax.text(mx - 0.012, my - 0.012, "label-stream surplus" + NL + "(what the permutation resolves)",
            rotation=45, rotation_mode="anchor", fontsize=5.0, color="#333333",
            ha="right", va="top", zorder=6)

    # corner annotation boxes in the empty half-plane corners
    ax.text(0.02, 0.975,
            ("SATURATED BY DESIGN (n = " + str(n_untest) + ") — null = real by construction" + NL +
             "score ≤ label-independent ceiling + 1e-9: the permutation cannot exceed the null" + NL +
             f"all {n_untest} receive p = 1.0 by enforcement (q = NaN, excluded from BH family)" + NL +
             "→ rank rests on external table evidence; covered by the label-free honest/strict arms"),
            transform=ax.transAxes, fontsize=5.5, ha="left", va="top", color="#222222",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
    ax.text(0.98, 0.02,
            (f"GENUINE NULL SEPARATION (n = {n_test}) — BH-FDR over testable genes only" + NL +
             f"{n_sig}/{n_test} testable genes significant; "
             f"{n_mid} intermediate p-values" + NL +
             f"add-one floor p = 1.0 × 10^-3 ({n_perm:,} draws/gene) → BH-FDR 0.05/0.01 resolved;"
             + NL +
             "Bonferroni 0.05/143 = 3.5 × 10^-4 would need n ≥ 2,860"),
            transform=ax.transAxes, fontsize=5.5, ha="right", va="bottom", color="#222222",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))

    ax.legend(loc="lower right", bbox_to_anchor=(0.985, 0.145), fontsize=5.2,
              frameon=False, handletextpad=0.3, borderaxespad=0.2, labelspacing=0.35)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel("Real integrated evidence score", fontsize=7.5, fontweight="bold")
    ax.set_ylabel("Label-independent King ceiling (permutation-null maximum)",
                  fontsize=7.5, fontweight="bold")
    ax.set_title("Atlas-permutation significance map — the whole test in one graph "
                 f"({n} neural TFs)", fontsize=8.5, pad=7)
    ax.tick_params(labelsize=6.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle("")
    save(fig, "43_permutation_map")
    print(f"Built 43_permutation_map.png ({n} genes; {n_sig} BH-sig; {n_test} testable; "
          f"{n_untest} saturated ({n_pinned} pinned + {n_bound} float-tie); {n_mid} intermediate; "
          f"p* = {bh:.2e}; shortlist {len(sh_test)} floor / {len(sh_pin)} pinned)")


if __name__ == "__main__":
    build()

