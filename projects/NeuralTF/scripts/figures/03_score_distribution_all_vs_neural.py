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

    fig, ax = plt.subplots(figsize=(W_15COL, 2.8))
    bins = np.linspace(0, max(a.max(), n.max()) * 1.02, 35)
    
    # Histograms
    ax.hist(a, bins=bins, density=True, alpha=0.35, color=C_ALL, edgecolor="none")
    ax.hist(n, bins=bins, density=True, alpha=0.55, color=C_A, edgecolor="none")

    # KDE overlay
    from scipy.stats import gaussian_kde
    x_grid = np.linspace(bins[0], bins[-1], 300)
    kde_a = gaussian_kde(a)
    kde_n = gaussian_kde(n)
    ax.plot(x_grid, kde_a(x_grid), color="#55606B", lw=1.2, label=f"All candidates (n = {len(a):,})")
    ax.plot(x_grid, kde_n(x_grid), color=C_A, lw=1.5, label=f"Neural TFs (n = {len(n):,})")

    ks, p = ks_2samp(a, n)
    p_str = "P < 10^{-30}" if p < 1e-30 else f"P = {p:.1e}"
    ax.text(0.96, 0.90, f"KS test: $D = {ks:.3f}$, ${p_str}$",
            transform=ax.transAxes, fontsize=6.2, ha="right", va="top", color="#333333")

    ax.set_xlabel("Integrated evidence score", fontsize=7.0)
    ax.set_ylabel("Probability density", fontsize=7.0)
    ax.set_title("Candidate score distribution (All TFs vs Neural-filtered)", fontsize=8.0, pad=6)

    ax.legend(frameon=False, fontsize=6.2, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save(fig, "03_score_distribution_all_vs_neural")

if __name__=="__main__": build()

