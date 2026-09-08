"""Score distribution: all vs neural-filtered candidates (histogram + KDE)."""
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

    fig, ax = plt.subplots(figsize=(W_15COL, 3.4))
    bins = np.linspace(0, max(a.max(), n.max()) * 1.02, 35)
    
    # Histograms
    ax.hist(a, bins=bins, density=True, alpha=0.35, color="#78909C", edgecolor="none")
    ax.hist(n, bins=bins, density=True, alpha=0.45, color=C_A, edgecolor="none")

    # KDE overlay
    from scipy.stats import gaussian_kde
    x_grid = np.linspace(bins[0], bins[-1], 300)
    kde_a = gaussian_kde(a)
    kde_n = gaussian_kde(n)
    ax.plot(x_grid, kde_a(x_grid), color="#455A64", lw=1.5, label=f"All candidates (n = {len(a):,})")
    ax.plot(x_grid, kde_n(x_grid), color=C_A, lw=1.8, label=f"Neural-enriched (n = {len(n):,})")

    ks, p = ks_2samp(a, n)
    p_str = "p < 10^{-30}" if p < 1e-30 else f"p = {p:.1e}"
    ax.text(0.96, 0.93, f"Kolmogorov–Smirnov\n$D = {ks:.3f}$\n${p_str}$",
            transform=ax.transAxes, fontsize=7.5, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.35", fc="#FAFAFA", ec="#CCCCCC", lw=0.6))

    ax.set_xlabel("Integrated evidence score", fontsize=8)
    ax.set_ylabel("Probability density", fontsize=8)
    ax.set_title("Neural filtering enriches for higher-scoring candidates",
                 fontweight="bold", fontsize=8.5, pad=8)

    ax.legend(frameon=False, fontsize=7.5, loc="upper center")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save(fig, "03_score_distribution_all_vs_neural")

if __name__=="__main__": build()

