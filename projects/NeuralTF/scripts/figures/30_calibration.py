"""Empirical decile calibration diagram for the integrated evidence score.

Single-panel evaluation of rank discrimination and score calibration:
- Validated neural TF rate across score deciles (D1 lowest to D10 highest)
- Wilson 95% binomial confidence intervals
- Genome-wide background prevalence reference line
- Top-decile enrichment fold-change annotation

The rate axis is kinked at 1% (see ``style.kinked_axis``): seven deciles sit at
exactly 0 and D8/D9 at 0.2 / 0.3%, all of which were 1-3 px slivers on the old
linear 0-8.6% scale set by the 5.3% top decile.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json

# kinked rate axis — [0, 1%] keeps full data resolution, the band
# above it is compressed; 1% sits above every background-level decile (0-0.35%)
# and above the prevalence line (0.58%) but well below the top decile (5.3%).
# KINK_FRAC = share of the axis height given to the linear segment below the kink.
KINK_PCT, KINK_FRAC = 1.0, 0.55


def build():
    data_path = RES / "calibration_stats.json"
    if not data_path.exists():
        raise FileNotFoundError(f"{data_path} missing — run scripts/stats/calibration.py first")

    with open(data_path) as f:
        data = json.load(f)

    stats_list = data.get("bin_stats", [])
    if not stats_list:
        raise ValueError("calibration_stats.json carries no bin_stats")
    bin_df = pd.DataFrame(stats_list)

    # Order deciles from D1 (lowest) to D10 (highest score)
    bin_df = bin_df.sort_values("mean_score", ascending=True).reset_index(drop=True)
    decile_labels = [f"D{i+1}\n({bin_df.iloc[i]['mean_score']:.2f})" for i in range(len(bin_df))]
    x = np.arange(len(bin_df))
    observed = bin_df["empirical_positive_rate"].to_numpy(dtype=float)
    counts = bin_df["n_candidates"].to_numpy(dtype=float)
    prevalence = float(data["prevalence"])

    fig, ax = plt.subplots(figsize=(W_15COL, 3.8), dpi=500)

    # Bars: Highlight bins enriched above background prevalence
    colors = [C_A if obs >= prevalence else "#D0D7DE" for obs in observed]
    bars = ax.bar(x, observed * 100, color=colors, width=0.62, edgecolor="none", zorder=3)

    # Wilson 95% confidence intervals per decile: the interval drawn is
    # [centre - half, centre + half] clipped to the valid 0-100% range,
    # not p_hat +/- half centred on the point estimate. A zero-rate decile
    # can only err upwards, so its CI is one-sided by construction.
    z = 1.96
    p_hat = np.clip(observed, 0.0, 1.0)
    denom = 1 + z**2 / counts
    centre = (p_hat + z**2 / (2 * counts)) / denom
    half = z * np.sqrt(p_hat * (1 - p_hat) / counts + z**2 / (4 * counts**2)) / denom
    ci_lo = np.clip(centre - half, 0.0, 1.0) * 100
    ci_hi = np.clip(centre + half, 0.0, 1.0) * 100
    ax.errorbar(x, observed * 100,
                yerr=np.vstack([observed * 100 - ci_lo, ci_hi - observed * 100]),
                fmt="none", ecolor="#333333", elinewidth=0.9, capsize=3, zorder=4,
                label="Wilson 95% CI")

    # Kinked rate axis: below 1% the axis keeps full resolution and
    # everything above is compressed, so the small deciles (D8/D9 at
    # 0.2/0.3% and every zero decile) do not shrink to a 1-3 px sliver the
    # way a linear 0-8.6% axis made them. A kink rather than a two-panel
    # break keeps every bar and error bar one intact shape. The callout
    # moves into the empty compressed band left of the top-decile bar.
    axis_hi = float(np.ceil(ci_hi[-1] * 1.10 * 10.0) / 10.0)
    annot_y = float(ci_hi[-1] * 0.90)

    # Prevalence reference line
    ax.axhline(y=prevalence * 100, color=C_HL, lw=1.1, linestyle="--", zorder=2,
                label=f"Genome prevalence ({prevalence*100:.2f}%)")

    # Top decile callout annotation with exact binomial test significance
    top_rate = observed[-1] * 100
    enrichment = observed[-1] / prevalence if prevalence > 0 else 0
    p_binom = data.get("top_decile_enrichment", {}).get("p_one_sided_binomial", None)
    if p_binom is not None:
        if p_binom < 1e-15:
            base, exp = f"{p_binom:.1e}".split("e")
            p_str = rf"P = {base} \times 10^{{{int(exp)}}}"
        else:
            p_str = f"P = {p_binom:.1e}"
        callout_txt = f"{top_rate:.1f}%\n({enrichment:.1f}\u00d7 enrichment,\n${p_str}$)"
    else:
        callout_txt = f"{top_rate:.1f}%\n({enrichment:.1f}\u00d7 enrichment)"

    ax.annotate(callout_txt,
                xy=(x[-1] - 0.33, top_rate),
                xytext=(x[-1] - 3.9, annot_y),
                arrowprops=dict(arrowstyle="->", color="#333333", lw=0.7),
                fontsize=6.0, ha="left", va="center", color="#222222")

    ax.set_xticks(x)
    ax.set_xticklabels(decile_labels, fontsize=6.2)
    ax.set_xlabel("Integrated Score Decile (Mean Score)", fontsize=7.0)
    ax.set_ylabel("RNAi-screened rate (%)", fontsize=7.0)
    kink_fac = kinked_axis(ax, "y", 0.0, KINK_PCT, axis_hi, frac_below=KINK_FRAC)
    yt = [t for t in (0, 0.5, KINK_PCT, 2, 3, 4, 5, 6, 7) if t <= axis_hi]
    ax.set_ylim(0, axis_hi)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{t:g}" for t in yt])
    # Disclose the non-linear axis at the kink itself (same wording as fig 06).
    ax.text(float(x[0]) - 0.42, KINK_PCT * 1.75,
            f"y-axis kinked at {KINK_PCT:g}%: {kink_fac:.1f}x compressed above",
            fontsize=5.4, color="#666666", ha="left", va="center")
    # header raised ~2.5 pt + subtitle pulled 3.5 pt closer to
    # the title; legend anchor dropped 0.864 -> 0.861 so the subtitle no
    # longer touches the legend row (was a ~2.3 pt bbox overlap).
    title_block(
        fig,
        "RNAi-Screened Enrichment Across Integrated-Score Deciles",
        f"Genome prevalence = {prevalence*100:.2f}% (dashed); top decile = {top_rate:.1f}% "
        f"({enrichment:.1f}x enriched)",
        y=0.9941, sub_y=0.9477,
    )
    from matplotlib.patches import Patch
    h30, l30 = ax.get_legend_handles_labels()
    h30.append(Patch(facecolor=C_A))
    l30.append("Decile ≥ prevalence")
    fig.legend(handles=h30, labels=l30, loc="lower center", bbox_to_anchor=(0.5, 0.861),
               frameon=False, fontsize=5.8, ncol=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.subplots_adjust(left=0.14, right=0.96, top=0.86, bottom=0.16)
    save(fig, "30_calibration")


if __name__ == "__main__":
    build()
