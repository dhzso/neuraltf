"""Integrated score distribution — all 249 TF candidates with top-10 overlaid."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np

def build():
    all249 = load_all()
    neural = load_neural()
    top10 = load_top10()

    fig, ax = plt.subplots(figsize=(7, 5))
    scores = all249["integrated_score"].dropna().values
    ax.hist(scores, bins=30, color=C_ALL, edgecolor="white", lw=0.3, alpha=0.6, label="All 249 TFs")
    n_scores = neural["integrated_score"].dropna().values
    ax.hist(n_scores, bins=25, color=C_NEURAL, edgecolor="white", lw=0.3, alpha=0.5, label="99 neural candidates")
    ax.axvline(np.median(scores), color="#999", ls="--", lw=0.8, label=f"Median (249) = {np.median(scores):.3f}")
    ax.axvline(np.median(n_scores), color=C_HL, ls=":", lw=0.8, label=f"Median (99) = {np.median(n_scores):.3f}")

    # Overlay top-10 positions
    top10_ids = set(top10["gene_id"].tolist())
    for track, c, sym in [("A", C_A, "D"), ("B", C_B, "s")]:
        tids = [g for g in top10["gene_id"] if g in top10_ids and top10[top10["gene_id"]==g].iloc[0].get("track","")==track]
        tdf = neural[neural["gene_id"].isin(tids)]
        if not tdf.empty:
            ax.scatter(tdf["integrated_score"], np.zeros(len(tdf))+0.5, c=c, marker=sym, s=50, zorder=5,
                       edgecolors="white", lw=0.5, label=f"Track {'A' if track=='A' else 'B'} top-10")
    ax.set_xlabel("Integrated evidence score")
    ax.set_ylabel("Count")
    ax.set_title("Integrated score distribution — 249 TF candidates", fontweight="bold", pad=8)
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "02_integrated_vs_composite")

if __name__=="__main__": build()
