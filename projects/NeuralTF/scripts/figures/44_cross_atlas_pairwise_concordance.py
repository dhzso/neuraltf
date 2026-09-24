"""Figure 44: Cross-Atlas Neural Effect Size Concordance — All Atlas Pairs.

Companion to figure 35. Figure 35 shows a single Fincher x Plass scatter; this
figure extends the exact same visualization to every pairwise atlas contrast
available in the pipeline's per-gene DE checkpoint (``de_pvalues.parquet``):

    panel a  Fincher et al. 2018  x  Plass et al. 2018
    panel b  Fincher et al. 2018  x  Cui et al. 2023
    panel c  Plass et al. 2018    x  Cui et al. 2023

All three atlases are Leiden-clustered and scored independently with the same
Wilcoxon best-cluster convention (``Pipeline.score_atlases`` ->
``_checkpoint_post_scoring``), so the per-gene
best-cluster ``log2`` fold changes are directly comparable across panels.

Conventions / data notes:
- Effect-size windows mirror figure 35 (Fincher -2..8, Plass -2..10); the Cui
  window (-2..10) retains >99% of its positive-lfc distribution. The pipeline
  records upregulation-only effect sizes, hence no negative fold changes.
- Uses the per-gene DE checkpoint only: a gene with effect sizes in both
  atlases of a pair already satisfies the >=2-atlas requirement of
  ``meta_analysis_pvalues.csv`` (verified: identical gene sets), so the meta
  table is not needed here.
- Neural candidates are the 143 genes of ``rank_neural.csv``
  (``neural_enriched > 0 OR rnai > 0``).
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

# (x-atlas, y-atlas, x label, y label)
PAIRS = [
    ("fincher", "plass", "Fincher et al.", "Plass et al."),
    ("fincher", "cui",   "Fincher et al.", "Cui et al."),
    ("plass",   "cui",   "Plass et al.",   "Cui et al."),
]

# Effect-size windows per atlas (identical to fig 35 for Fincher/Plass; the
# Cui window retains >99% of its positive-lfc mass).
LFC_WINDOW = {"fincher": (-2.0, 8.0), "plass": (-2.0, 10.0), "cui": (-2.0, 10.0)}
# Axis limits per atlas (mirrors fig 35's asymmetric windows).
AXIS_LIMIT = {"fincher": (-1.5, 9.5), "plass": (-1.5, 10.5), "cui": (-1.5, 10.5)}

# Curated per-panel labels: every plotted value is read from the data; the
# offsets are layout only, and a gene absent from a panel is skipped (the
# same contract as fig 35's annotation block).
ANNOTATIONS = {
    ("fincher", "plass"): [
        ("dd_Smed_v6_1854_0_1",  "dd1854",  (-18, 12), "right",  "bottom"),
        ("dd_Smed_v6_10038_0_1", "Zeb-1",   (-14, 14), "right",  "bottom"),
        ("dd_Smed_v6_11150_0_1", "dd11150", (-18, 2),  "right",  "center"),
    ],
    ("fincher", "cui"): [
        ("dd_Smed_v6_14753_0_1", "ascl-2",  (0, 12),   "center", "bottom"),
        ("dd_Smed_v6_1854_0_1",  "dd1854",  (-16, 12), "right",  "bottom"),
        ("dd_Smed_v6_17726_0_1", "Pax6A",   (-14, 8),  "right",  "bottom"),
    ],
    ("plass", "cui"): [
        ("dd_Smed_v6_17143_0_1", "tbx2/3b", (12, 8),   "left",   "bottom"),
        ("dd_Smed_v6_10038_0_1", "Zeb-1",   (12, 6),   "left",   "bottom"),
        ("dd_Smed_v6_2946_0_1",  "dd2946",  (-16, 4),  "right",  "center"),
    ],
}



def build():
    de_path = RUN / "de_pvalues.parquet"
    if not de_path.exists():
        raise FileNotFoundError(
            f"{de_path} missing — re-run the pipeline (scripts/run.py); "
            "checkpoint 03 writes the per-atlas DE p/lfc columns."
        )
    de = pd.read_parquet(de_path)
    neural = load_neural()

    fig, axes = plt.subplots(1, 3, figsize=(W_2COL, 3.6), dpi=500)
    per_panel = []  # (label, n, spearman, pearson) — feeds the subtitle

    for ax, letter, (ax_key, ay_key, ax_lab, ay_lab) in zip(axes, "abc", PAIRS):
        x_col, y_col = f"{ax_key}_lfc", f"{ay_key}_lfc"
        sub = de.dropna(subset=[x_col, y_col]).copy()
        sub = sub[sub[x_col].between(*LFC_WINDOW[ax_key])
                  & sub[y_col].between(*LFC_WINDOW[ay_key])]
        if len(sub) < 100:
            raise ValueError(
                f"{ax_key} x {ay_key}: only {len(sub)} genes with usable "
                "effect sizes — regenerate de_pvalues.parquet."
            )

        rho, _ = spearmanr(sub[x_col], sub[y_col])
        r_p, _ = pearsonr(sub[x_col], sub[y_col])
        slope, intercept = np.polyfit(sub[x_col], sub[y_col], 1)
        name = f"{ax_lab.split()[0]} × {ay_lab.split()[0]}"
        per_panel.append((name, len(sub), rho, r_p))

        # Reference lines: zero dotted, identity dashed, least-squares fit
        ax.axhline(0, color="#D0D7DE", lw=0.6, ls=":", zorder=1)
        ax.axvline(0, color="#D0D7DE", lw=0.6, ls=":", zorder=1)
        ax.plot(AXIS_LIMIT[ax_key], AXIS_LIMIT[ay_key], color="#888888",
                lw=0.8, ls="--", zorder=2, label="Identity ($y = x$)")

        x_vals = np.linspace(sub[x_col].min(), sub[x_col].max(), 100)
        ax.plot(x_vals, slope * x_vals + intercept, color=C_A, lw=1.0, zorder=3,
                label="Linear fit")

        ax.scatter(sub[x_col], sub[y_col], c=C_A, s=5, alpha=0.22,
                   edgecolor="none", rasterized=True, zorder=3,
                   label="All genes (both atlases)")

        neural_common = neural.merge(sub[["v6_id", x_col, y_col]],
                                     left_on="gene_id", right_on="v6_id",
                                     how="inner")
        ax.scatter(neural_common[x_col], neural_common[y_col], c=C_B,
                   edgecolor="#222222", linewidth=0.6, s=26, zorder=6,
                   label="Neural TF candidates")

        for gid, nm, offset, ha, va in ANNOTATIONS.get((ax_key, ay_key), []):
            row = neural_common[neural_common["v6_id"] == gid]
            if len(row) == 0:
                continue
            ax.annotate(
                nm,
                xy=(row.iloc[0][x_col], row.iloc[0][y_col]),
                xytext=offset,
                textcoords="offset points",
                fontsize=5.6,
                color="#111111",
                ha=ha,
                va=va,
                arrowprops=dict(arrowstyle="-", color="#555555", lw=0.5),
                zorder=7,
            )

        ax.set_xlim(*AXIS_LIMIT[ax_key])
        ax.set_ylim(*AXIS_LIMIT[ay_key])
        ax.set_xlabel(f"{ax_lab} $\\log_2$ fold change", fontsize=6.4)
        ax.set_ylabel(f"{ay_lab} $\\log_2$ fold change", fontsize=6.4)
        ax.set_title(f"{name}\n$n$ = {len(sub):,}, $\\rho$ = {rho:.2f}, $r$ = {r_p:.2f}",
                     fontsize=7.2, pad=5)
        panel_tag(ax, letter, x=-0.17, y=1.09, fontsize=8.0)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.850),
               ncol=4, frameon=False, fontsize=6.0, handletextpad=0.4,
               columnspacing=1.4)

    # Subtitle numbers are computed from the data (never hard-coded).
    sub_txt = "; ".join(f"{nm} $\\rho$ = {rho:.2f}" for nm, _, rho, _ in per_panel)
    title_block(
        fig,
        "Cross-Atlas Neural Effect Size Concordance — All Atlas Pairs",
        "Best-cluster $\\log_2$ fold changes (Spearman): " + sub_txt,
        y=0.995, sub_y=0.9325,
    )
    fig.subplots_adjust(left=0.085, right=0.985, top=0.775, bottom=0.19,
                        wspace=0.22)
    save(fig, "44_cross_atlas_pairwise_concordance")


if __name__ == "__main__":
    build()
