"""Integrated score distribution with composite scores for top-10 overlaid."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd

def build():
    neural = load_neural()
    top10 = load_top10()
    centered = load_centered()
    uniform = load_uniform()

    # Build short gene IDs
    def short(gid):
        parts = gid.split("_")
        for p in parts:
            if p.startswith("dd") and p[2:].isdigit(): return p
        return gid

    fig, ax = plt.subplots(figsize=(7, 5))
    scores = neural["integrated_score"].values
    ax.hist(scores, bins=25, color=C_ALL, edgecolor="white", lw=0.3, alpha=0.6, label="All 99 neural")
    ax.axvline(np.median(scores), color="#999", ls="--", lw=0.8, label=f"Median = {np.median(scores):.3f}")

    # Overlay top-10: show integrated_score from rank_neural AND composite from top10
    top10_ids = set(top10["gene_id"].tolist())
    for track, c, sym in [("A", C_A, "D"), ("B", C_B, "s")]:
        tids = [g for g in top10["gene_id"] if g in top10_ids and top10[top10["gene_id"]==g].iloc[0].get("track","")==track]
        tdf = neural[neural["gene_id"].isin(tids)]
        if not tdf.empty:
            ax.scatter(tdf["integrated_score"], np.zeros(len(tdf))+0.3, c=c, marker=sym, s=40, zorder=5,
                       edgecolors="white", lw=0.5, label=f"Track {'A' if track=='A' else 'B'} (integrated)")
    ax.set_xlabel("Integrated evidence score"); ax.set_ylabel("Count")
    ax.set_title("Integrated score distribution — 99 neural candidates", fontweight="bold", pad=8)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "02_integrated_vs_composite")

if __name__=="__main__": build()
