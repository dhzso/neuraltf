"""Perez influence comparison — TF lineage class vs neuron influence scores.

Reads the Perez tables from the pipeline run:
  Panel a: Per-gene integrated scores grouped by Perez TF lineage class
           (neural-class, other-class, unclassified).
  Panel b: Distribution of the Perez ANANSE neuron influence stream
           for RBH-mapped candidates.
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    rank = load_all()

    has_lineage = "perez_lineage" in rank.columns and rank["perez_lineage"].notna().any()
    has_infl = "perez_influence" in rank.columns and rank["perez_influence"].notna().any()
    if not (has_lineage or has_infl):
        raise FileNotFoundError(
            "rank.csv carries no perez_lineage/perez_influence values."
        )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_2COL, 2.7))

    # Panel a: neural-class vs other-class TFs (perez_lineage stream)
    if has_lineage:
        neural_cls = rank[rank["perez_lineage"] == 1.0]["integrated_score"].dropna()
        other_cls = rank[rank["perez_lineage"] == 0.5]["integrated_score"].dropna()
        absent = rank[rank["perez_lineage"] == 0.0]["integrated_score"].dropna()
        data, labels, colors = [], [], []
        for vals, lab, col in (
            (neural_cls, "Neural-class TFs", C_A),
            (other_cls, "Other-class TFs", C_B),
            (absent, "Unclassified", "#CCCCCC"),
        ):
            if len(vals) > 0:
                data.append(vals.values)
                labels.append(f"{lab}\n(n={len(vals):,})")
                colors.append(col)
        if data:
            bp = ax1.boxplot(data, patch_artist=True, widths=0.55,
                             medianprops=dict(color="#111111", lw=1.2),
                             boxprops=dict(lw=0.7),
                             whiskerprops=dict(lw=0.7, color="#555555"),
                             capprops=dict(lw=0.7, color="#555555"),
                             flierprops=dict(marker=".", markersize=2, alpha=0.3))
            for patch, col in zip(bp["boxes"], colors):
                patch.set_facecolor(col)
                patch.set_alpha(0.7)
                patch.set_edgecolor("#333333")
            ax1.set_xticklabels(labels, fontsize=7)
        ax1.set_ylabel("Integrated score", fontsize=7.5)
        ax1.set_title("Score by lineage classification", fontsize=8, pad=4)
    else:
        ax1.text(0.5, 0.5, "perez_lineage stream empty in this run",
                 ha="center", va="center", transform=ax1.transAxes,
                 fontsize=8, color="#999999")
    panel_tag(ax1, "a")

    # Panel b: perez_influence stream distribution
    if has_infl:
        infl = rank["perez_influence"].dropna()
        infl_nz = infl[infl > 0].values
        ax2.hist(infl_nz, bins=25, color=STREAM_C["perez_influence"],
                 alpha=0.75, edgecolor="none")
        med = np.median(infl_nz) if len(infl_nz) > 0 else 0
        ax2.axvline(x=med, color=C_HL, lw=1.2, linestyle="--",
                    label=f"Median ({med:.2f})")
        ax2.set_xlabel("ANANSE neuron influence score", fontsize=7.5)
        ax2.set_ylabel("Number of candidates", fontsize=7.5)
        ax2.set_title("Influence score distribution", fontsize=8, pad=4)
        ax2.legend(loc="upper right", frameon=False, fontsize=7)
    else:
        ax2.text(0.5, 0.5, "perez_influence stream empty in this run",
                 ha="center", va="center", transform=ax2.transAxes,
                 fontsize=8, color="#999999")
    panel_tag(ax2, "b")

    fig.tight_layout()
    save(fig, "32_perez_influence_comparison")

if __name__ == "__main__":
    build()

