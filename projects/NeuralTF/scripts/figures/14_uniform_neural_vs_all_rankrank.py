"""Fixed-weight vs Uniform Dirichlet — all TF candidates (rank-rank comparison)."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np
from scipy.stats import spearmanr

def build():
    all_cand = load_all()
    unif = load_uniform_full()
    top10 = load_uniform()
    top10_ids = set(top10["gene_id"].tolist())
    track_map = dict(zip(top10["gene_id"], top10.get("track", [""]*len(top10))))

    all_cand["r_fixed"] = all_cand["integrated_score"].rank(ascending=False)
    unif["r_unif"] = unif["uniform_median_score"].rank(ascending=False)
    merged = all_cand[["gene_id", "r_fixed"]].merge(unif[["gene_id", "r_unif"]], on="gene_id", how="inner")

    fig, ax = plt.subplots(figsize=(W_15COL, 3.8))
    non_top = merged[~merged["gene_id"].isin(top10_ids)]
    ax.scatter(non_top["r_fixed"], non_top["r_unif"], s=8, c="#B0BEC5", alpha=0.3, edgecolors="none", label="Other candidates")
    for track, c in [("A", C_A), ("B", C_B)]:
        tids = [g for g, t in track_map.items() if t == track]
        sub = merged[merged["gene_id"].isin(tids)]
        if not sub.empty:
            ax.scatter(sub["r_fixed"], sub["r_unif"], s=32, c=c, edgecolors="white", lw=0.6, zorder=5, label=f"Track {track}")
            
    # Add clean leader annotations for top candidates
    for _, r in merged[merged["gene_id"].isin(top10_ids)].iterrows():
        nm = label(top10, r["gene_id"])
        ax.annotate(nm, (r["r_fixed"], r["r_unif"]), fontsize=6, ha="left", va="bottom",
                    xytext=(4, 2), textcoords="offset points", fontweight="bold",
                    color=C_A if track_map.get(r["gene_id"]) == "A" else C_B)
                    
    lim = max(merged["r_fixed"].max(), merged["r_unif"].max()) * 1.05
    ax.plot([0, lim], [0, lim], "--", color="#555555", lw=0.8, label="y = x (identity)")
    rho, p = spearmanr(merged["r_fixed"], merged["r_unif"])
    p_str = "p < 10^{-30}" if p < 1e-30 else f"p = {p:.1e}"
    ax.text(0.95, 0.08, f"$r_s = {rho:.3f}$\n${p_str}$\n(n = {len(merged):,})",
            transform=ax.transAxes, fontsize=7.5, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", fc="#FAFAFA", ec="#CCCCCC", lw=0.6))
            
    ax.set_xlabel(f"Fixed-weight integrated rank (n = {len(merged):,})", fontsize=8)
    ax.set_ylabel("Uniform Dirichlet median rank (α = 1)", fontsize=8)
    ax.set_title("Rank conservation across candidate universe under non-informative Dirichlet weights",
                 fontweight="bold", fontsize=8.5, pad=8)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "14_uniform_neural_vs_all_rankrank")

if __name__ == "__main__":
    build()

