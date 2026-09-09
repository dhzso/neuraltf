"""Top 10 prioritized candidates — unified single-panel horizontal atlas."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# TF family colors (from InterPro domain annotations)
TF_FAMILIES = {
    "Homeobox": "#0072B2",
    "bHLH": "#E69F00",
    "Zinc finger": "#009E73",
    "Fork head": "#D55E00",
    "T-box": "#CC79A7",
    "p53": "#56B4E9",
}


def tf_family_color(domains_str):
    d = str(domains_str)
    for fam, c in TF_FAMILIES.items():
        if fam.lower() in d.lower():
            return fam, c
    return "Other", "#AAAAAA"


def build():
    neural = load_neural()
    top10 = load_top10()

    if top10.empty:
        fig, ax = plt.subplots(figsize=(W_15COL, 4.0))
        ax.text(0.5, 0.5, "No candidates", ha="center", va="center")
        save(fig, "05_top10_candidate_atlas")
        return

    # Extract candidate metadata
    records = []
    for _, row in top10.iterrows():
        gid = row["gene_id"]
        track = row.get("track", "B")
        nm = label(neural, gid)
        n_row = neural[neural["gene_id"] == gid]
        base = n_row.iloc[0].get("integrated_score", np.nan) if len(n_row) > 0 else np.nan
        comp = row.get("composite_score", np.nan)
        orth = str(row.get("human_ortholog", "") or "")
        orth_clean = orth if orth.lower() not in ("nan", "none", "") else "—"
        domains = str(row.get("interpro_domains", row.get("domains_all", "")))
        fam, _ = tf_family_color(domains)
        records.append({
            "gene_id": gid,
            "name": nm,
            "track": track,
            "base": base,
            "composite": comp,
            "ortholog": orth_clean,
            "family": fam,
        })
    df = pd.DataFrame(records)

    # Sort Track A and Track B individually by composite score ascending
    sub_a = df[df["track"] == "A"].sort_values("composite", ascending=True)
    sub_b = df[df["track"] == "B"].sort_values("composite", ascending=True)

    # Combined layout: Track B on bottom (indices 0..4), Track A on top (indices 6..10), gap at 5
    n_b = len(sub_b)
    n_a = len(sub_a)
    y_b = np.arange(n_b)
    gap = 1.2
    y_a = n_b + gap + np.arange(n_a)

    fig, ax = plt.subplots(figsize=(W_15COL, 4.6), dpi=500)

    # Bars for Track B
    ax.barh(y_b, sub_b["composite"], height=0.62, color=C_B, alpha=0.88,
            edgecolor="none", label="Track B (Candidate)")
    ax.scatter(sub_b["base"], y_b, color="#222222", s=28, zorder=4,
               edgecolors="white", lw=0.6, label="Base score")

    for i, (_, r) in enumerate(sub_b.iterrows()):
        ax.text(r["composite"] + 0.015, y_b[i],
                f"{r['composite']:.2f} ({r['family']})",
                fontsize=6.8, va="center", color="#222222")

    # Bars for Track A
    ax.barh(y_a, sub_a["composite"], height=0.62, color=C_A, alpha=0.88,
            edgecolor="none", label="Track A (Benchmark)")
    ax.scatter(sub_a["base"], y_a, color="#222222", s=28, zorder=4,
               edgecolors="white", lw=0.6)

    for i, (_, r) in enumerate(sub_a.iterrows()):
        ax.text(r["composite"] + 0.015, y_a[i],
                f"{r['composite']:.2f} ({r['family']})",
                fontsize=6.8, va="center", color="#222222")

    # Divider line and section labels
    div_y = n_b + gap / 2.0 - 0.5
    ax.axhline(div_y, color="#CCCCCC", lw=0.8, ls="--")

    # Y-ticks
    all_y = np.concatenate([y_b, y_a])
    all_names = list(sub_b["name"]) + list(sub_a["name"])
    ax.set_yticks(all_y)
    ax.set_yticklabels(all_names, fontsize=7.5)

    # Section annotations on the right or near divider
    ax.text(0.01, div_y + 0.25, "TRACK A: RNAi-Validated Regulators",
            fontsize=7.0, fontweight="bold", color=C_A, va="bottom")
    ax.text(0.01, div_y - 0.25, "TRACK B: Novel Candidate Regulators",
            fontsize=7.0, fontweight="bold", color=C_B, va="top")

    ax.set_xlabel("Prioritization Score", fontsize=8, fontweight="bold")
    ax.set_xlim(0, 1.48)
    ax.set_ylim(-0.8, y_a[-1] + 0.8)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Clean unified legend placed cleanly above axes
    leg_handles = [
        Patch(facecolor=C_A, label="Track A (RNAi-validated)"),
        Patch(facecolor=C_B, label="Track B (Novel candidate)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#222222",
               markeredgecolor="white", markersize=5, label="Base score"),
    ]
    ax.legend(handles=leg_handles, loc="lower center", bbox_to_anchor=(0.5, 1.02),
              ncol=3, frameon=False, fontsize=7)

    fig.suptitle("Prioritized Neural Transcription Factors (Top 10 Candidates)",
                 fontweight="bold", fontsize=8.5, y=0.98)
    fig.subplots_adjust(left=0.18, right=0.96, top=0.88, bottom=0.12)
    save(fig, "05_top10_candidate_atlas")


if __name__ == "__main__":
    build()
