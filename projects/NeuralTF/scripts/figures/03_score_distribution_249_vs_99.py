"""Score distribution: 249 full vs 99 neural-filtered candidates (histogram + KDE)."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np
from scipy.stats import ks_2samp

def build():
    all_df = load_all()
    neural = load_neural()
    col = "integrated_score"
    a = all_df[col].dropna().values
    n = neural[col].dropna().values

    fig, ax = plt.subplots(figsize=(7, 4.5))
    bins = np.linspace(0, max(a.max(),n.max())*1.02, 35)
    ax.hist(a, bins=bins, density=True, alpha=0.5, color=C_ALL, edgecolor="white", lw=0.3, label=f"All candidates (n={len(a)})")
    ax.hist(n, bins=bins, density=True, alpha=0.6, color=C_CENTERED, edgecolor="white", lw=0.3, label=f"Neural-enriched (n={len(n)})")

    # KDE overlay
    from scipy.stats import gaussian_kde
    x_grid = np.linspace(bins[0], bins[-1], 200)
    for vals, c, lbl in [(a, C_ALL, f"All candidates KDE (n={len(a)})"), (n, C_CENTERED, f"Neural KDE (n={len(n)})")]:
        kde = gaussian_kde(vals)
        ax.plot(x_grid, kde(x_grid), color=c, lw=1.2, label=lbl)

    ks, p = ks_2samp(a, n)
    ax.text(0.97, 0.95, f"KS = {ks:.3f}\np = {p:.1e}", transform=ax.transAxes,
            fontsize=8, ha="right", va="top", bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
    ax.set_xlabel("Integrated evidence score"); ax.set_ylabel("Relative frequency of candidates")
    ax.set_title("Neural filtering enriches for higher-scoring candidates\n"
                 "(Neural-enriched vs All candidates; KS test shown)", fontweight="bold", pad=8)

    ax.legend(frameon=False, fontsize=7, loc="upper left"); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "03_score_distribution_249_vs_99")

if __name__=="__main__": build()
