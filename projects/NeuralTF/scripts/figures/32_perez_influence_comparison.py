"""Perez influence comparison — TF lineage class vs neuron influence scores.

Reads the REAL Perez tables from the pipeline run (no invented columns):
  Panel A: per-gene perez_lineage stream scores grouped by neural vs
           other TF class (from rank.csv, produced by integrate_perez).
  Panel B: distribution of the perez_influence stream (MOESM19 neuron
           sheet, RBH-mapped), for genes that have one.
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    rank = load_all()  # raises if missing — no silent fallback

    # need at least one of the two Perez streams
    has_lineage = "perez_lineage" in rank.columns and rank["perez_lineage"].notna().any()
    has_infl = "perez_influence" in rank.columns and rank["perez_influence"].notna().any()
    if not (has_lineage or has_infl):
        raise FileNotFoundError(
            "rank.csv carries no perez_lineage/perez_influence values — "
            "re-run the pipeline (the current run predates the Perez streams)."
        )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    # Panel A: neural-class vs other-class TFs (perez_lineage stream)
    if has_lineage:
        neural_cls = rank[rank["perez_lineage"] == 1.0]["integrated_score"].dropna()
        other_cls = rank[rank["perez_lineage"] == 0.5]["integrated_score"].dropna()
        absent = rank[rank["perez_lineage"] == 0.0]["integrated_score"].dropna()
        data, labels, colors = [], [], []
        for vals, lab, col in (
            (neural_cls, "Neural-class TFs", C_A),
            (other_cls, "Other-class TFs", C_B),
            (absent, "Not in Perez", "#BBBBBB"),
        ):
            if len(vals) > 0:
                data.append(vals.values)
                labels.append(f"{lab}\n(n={len(vals)})")
                colors.append(col)
        if data:
            bp = ax1.boxplot(data, patch_artist=True,
                             medianprops=dict(color=C_HL))
            for patch, col in zip(bp["boxes"], colors):
                patch.set_facecolor(col)
                patch.set_alpha(0.6)
            ax1.set_xticklabels(labels, fontsize=7)
        ax1.set_ylabel("Integrated score")
        ax1.set_title("Integrated score by Perez TF lineage class",
                      fontweight="bold")
    else:
        ax1.text(0.5, 0.5, "perez_lineage stream empty in this run",
                 ha="center", va="center", transform=ax1.transAxes,
                 fontsize=9, color="#999999")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: perez_influence stream distribution (MOESM19 neuron sheet)
    if has_infl:
        infl = rank["perez_influence"].dropna().values
        ax2.hist(infl, bins=25, color=STREAM_C["perez_influence"],
                 alpha=0.65, edgecolor="white")
        ax2.axvline(x=np.median(infl), color=C_HL, lw=1.5, linestyle="--",
                    label=f"Median = {np.median(infl):.3f}")
        ax2.set_xlabel("Perez ANANSE neuron influence score")
        ax2.set_ylabel("Number of candidates")
        ax2.set_title(f"Neuron-fate regulatory influence\n({len(infl)} candidates with RBH-mapped influence)",
                      fontweight="bold")
        ax2.legend(fontsize=7)
    else:
        ax2.text(0.5, 0.5, "perez_influence stream empty in this run",
                 ha="center", va="center", transform=ax2.transAxes,
                 fontsize=9, color="#999999")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "32_perez_influence_comparison")

if __name__ == "__main__":
    build()
