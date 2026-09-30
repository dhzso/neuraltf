"""Figure 39: Score-Uncertainty Forest Plot.

Horizontal bands showing the 95% weight-perturbation range of each top-10
candidate under Dirichlet weight uncertainty, for BOTH priors the pipeline
computes: centered (k=40) and uniform (alpha=1). Candidates are ordered by
integrated score (highest at top).

Both priors are drawn, because scripts/stats/bootstrap_confidence.py
writes both arms into results/bootstrap_scores_ci.csv (centered_* AND
uniform_* columns) and the uniform prior — the one the rank-distribution
analysis in fig 06 uses — has bands ~2x wider on the shortlist
(0.18-0.50 vs 0.09-0.21). Each candidate carries two bands: centered
solid on top, uniform dashed below, joined by a tick at the fixed-weight
integrated score.

Each band also carries its draw median (hollow marker — circle =
centered arm, square = uniform arm). The band alone shows only spread,
and the vertical fixed-weight tick is a different quantity from the
centre of the draw distribution. Marker values come from the
centered_median / uniform_median columns written by
scripts/stats/bootstrap_confidence.py; an older CSV falls back to the
corresponding *_mean column.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# (column prefix, row offset in row units, line style, line width, median marker)
ARMS = (("centered", +0.20, "-", 1.8, "o"),
        ("uniform", -0.20, (0, (3.2, 1.7)), 1.4, "s"))


def build():
    ci = pd.read_csv(RES / "bootstrap_scores_ci.csv")
    top10 = load_top10()

    # Merge to get track info and limit to top-10
    top10_ids = set(top10["gene_id"].values)
    ci_top = ci[ci["gene_id"].isin(top10_ids)].copy()
    ci_top = ci_top.merge(top10[["gene_id", "track"]], on="gene_id", how="left")
    ci_top = ci_top.sort_values("integrated_score", ascending=True)  # bottom-to-top

    # Build labels
    labels = []
    for _, row in ci_top.iterrows():
        name = label(ci_top, row["gene_id"])
        track = row.get("track", "?")
        labels.append(f"{name}  (Track {track})")

    # Axis range follows BOTH arms: the uniform bands reach far below the
    # centered ones (0.46 vs 0.69 on the shortlist).
    band_lo = float(ci_top[[f"{k}_ci_95_lo" for k, *_ in ARMS]].min().min())
    band_hi = float(ci_top[[f"{k}_ci_95_hi" for k, *_ in ARMS]].max().max())
    x_lo = band_lo - 0.03
    x_hi = band_hi + 0.10  # room for the two per-band Delta annotations

    fig, ax = plt.subplots(figsize=(W_1COL, 4.0), dpi=500)

    # Color by track
    track_colors = {"A": C_A, "B": C_HL}

    for i, row in enumerate(ci_top.itertuples()):
        c = track_colors.get(getattr(row, "track", None), "#888888")
        # Fixed-weight score: one tick spanning both of this gene's bands
        # (identical under either prior; it sits inside both bands for all 10).
        ax.plot([row.integrated_score, row.integrated_score], [i - 0.30, i + 0.30],
                color="#333333", lw=0.8, zorder=2)
        for key, dy, ls, lw, mk in ARMS:
            lo = getattr(row, f"{key}_ci_95_lo")
            hi = getattr(row, f"{key}_ci_95_hi")
            ax.plot([lo, hi], [i + dy, i + dy], color=c, lw=lw, ls=ls,
                    solid_capstyle="round", zorder=3)
            # Draw median: the band's location (falls back to mean if absent)
            med = getattr(row, f"{key}_median", None)
            if med is None:
                med = getattr(row, f"{key}_mean", None)
            if med is not None and np.isfinite(med):
                ax.plot([med], [i + dy], marker=mk, ms=3.2, mfc="white",
                        mec=c, mew=0.9, ls="none", zorder=4)
            ax.text(hi + 0.008, i + dy, f"\u0394 {hi - lo:.3f}", va="center",
                    fontsize=5.2, color="#555555")

    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlabel("Integrated evidence score", fontsize=7,
                  fontweight="bold")
    ax.set_xlim(x_lo, x_hi)
    ax.set_xticks([t for t in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0) if x_lo < t < x_hi])
    ax.set_ylim(-0.55, len(labels) - 0.45)
    sub_bot = title_block(
        fig,
        "Weight-Perturbation Bands (Dirichlet Draws)",
        "Bands = 2.5-97.5% range over 1,000 Dirichlet weight draws; not CIs",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend above the axes, anchored to the subtitle so it stays centred on the
    # FIGURE (an axes-level legend is centred on the axes and, at 5 entries,
    # overhung the 3.5 in canvas by 0.12 in). The legend fills column-major, so
    # the order puts the two track colours in column 1, the two arm styles in
    # column 2 and the point reference in column 3.
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=C_A, lw=2, label="Track A (RNAi-screened)"),
        Line2D([0], [0], color=C_HL, lw=2, label="Track B (not tested)"),
        Line2D([0], [0], color="#555555", lw=1.8, label="Centered k=40 (solid)"),
        Line2D([0], [0], color="#555555", lw=1.4, ls=(0, (3.2, 1.7)),
               label="Uniform \u03b1=1 (dashed)"),
        Line2D([0], [0], marker="|", color="#333333", lw=0, markersize=7,
               markeredgewidth=0.8, label="Fixed-weight score"),
        Line2D([0], [0], color="none", marker="o", markersize=3.2,
               mfc="white", mec="#555555", mew=0.9, label="Draw median"),
    ]
    _pt = 1.0 / (72.0 * fig.get_size_inches()[1])
    leg_top = sub_bot - 6 * _pt  # 6 pt gap below the subtitle
    fig.legend(handles=legend_elements, frameon=False, fontsize=5.2, loc="upper center",
               bbox_to_anchor=(0.5, leg_top), ncol=3)

    fig.subplots_adjust(left=0.32, right=0.92, top=0.80, bottom=0.12)
    save(fig, "39_score_uncertainty_forest")
    print("Built 39_score_uncertainty_forest.png")


if __name__ == "__main__":
    build()
