"""Uniform Dirichlet — 99 neural vs 249 wide (rank-rank comparison)."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np
from scipy.stats import spearmanr

def build():
    u99 = load_unif99()
    u249 = load_unif249()
    top10 = load_top10()
    top10_ids = set(top10["gene_id"].tolist())
    track_map = dict(zip(top10["gene_id"], top10.get("track",[""]*len(top10))))

    u99["r99"] = u99["uniform_median_score"].rank(ascending=False)
    u249["r249"] = u249["uniform_median_score"].rank(ascending=False)
    merged = u99[["gene_id","r99"]].merge(u249[["gene_id","r249"]], on="gene_id", how="inner")

    fig, ax = plt.subplots(figsize=(6, 5.5))
    non_top = merged[~merged["gene_id"].isin(top10_ids)]
    ax.scatter(non_top["r99"], non_top["r249"], s=10, c=C_NEURAL, alpha=0.3, edgecolors="none", label="Other neural")
    for track, c in [("A",C_A),("B",C_B)]:
        tids = [g for g,t in track_map.items() if t==track]
        sub = merged[merged["gene_id"].isin(tids)]
        if not sub.empty:
            ax.scatter(sub["r99"], sub["r249"], s=40, c=c, edgecolors="white", lw=0.5, zorder=5, label=f"Track {track}")
    for _, r in merged[merged["gene_id"].isin(top10_ids)].iterrows():
        ax.annotate(label(top10, r["gene_id"]), (r["r99"], r["r249"]), fontsize=6, ha="center", va="bottom",
                    xytext=(0,4), textcoords="offset points", fontweight="bold")
    lim = max(merged["r99"].max(), merged["r249"].max())*1.05
    ax.plot([0,lim],[0,lim],"--",color="#999",lw=0.8)
    rho, p = spearmanr(merged["r99"], merged["r249"])
    ax.text(0.95, 0.05, f"Spearman rho = {rho:.3f}\np = {p:.1e}", transform=ax.transAxes,
            fontsize=8, ha="right", va="bottom", bbox=dict(boxstyle="round",fc="white",alpha=0.8))
    ax.set_xlabel(f"Rank ({len(u99)} neural candidates, uniform Dirichlet)")
    ax.set_ylabel(f"Rank ({len(u249)} full candidates, uniform Dirichlet)")
    ax.set_title(f"Neural-filtered rankings recover full-universe rankings (rho = {rho:.2f})",
                 fontweight="bold", pad=8)
    ax.legend(frameon=False, fontsize=7, loc="upper left"); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    fig.tight_layout(); save(fig, "14_uniform_99vs249_rankrank")

if __name__=="__main__": build()
