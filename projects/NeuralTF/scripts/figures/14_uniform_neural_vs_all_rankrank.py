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

    fig, ax = plt.subplots(figsize=(6, 5.5))
    non_top = merged[~merged["gene_id"].isin(top10_ids)]
    ax.scatter(non_top["r_fixed"], non_top["r_unif"], s=10, c=C_NEURAL, alpha=0.3, edgecolors="none", label="Other candidates")
    for track, c in [("A", C_A), ("B", C_B)]:
        tids = [g for g, t in track_map.items() if t == track]
        sub = merged[merged["gene_id"].isin(tids)]
        if not sub.empty:
            ax.scatter(sub["r_fixed"], sub["r_unif"], s=40, c=c, edgecolors="white", lw=0.5, zorder=5, label=f"Track {track}")
    for _, r in merged[merged["gene_id"].isin(top10_ids)].iterrows():
        ax.annotate(label(top10, r["gene_id"]), (r["r_fixed"], r["r_unif"]), fontsize=6, ha="center", va="bottom",
                    xytext=(0, 4), textcoords="offset points", fontweight="bold")
    lim = max(merged["r_fixed"].max(), merged["r_unif"].max()) * 1.05
    ax.plot([0, lim], [0, lim], "--", color="#999", lw=0.8)
    rho, p = spearmanr(merged["r_fixed"], merged["r_unif"])
    ax.text(0.95, 0.05, f"Spearman rho = {rho:.3f}\np = {p:.1e}", transform=ax.transAxes,
            fontsize=8, ha="right", va="bottom", bbox=dict(boxstyle="round", fc="white", alpha=0.8))
    ax.set_xlabel(f"Fixed-weight integrated rank (n = {len(merged)})")
    ax.set_ylabel("Uniform Dirichlet rank (non-informative prior)")
    ax.set_title(f"Rank stability under uniform Dirichlet perturbation (rho = {rho:.2f})",
                 fontweight="bold", pad=8)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "14_uniform_neural_vs_all_rankrank")

if __name__ == "__main__":
    build()
