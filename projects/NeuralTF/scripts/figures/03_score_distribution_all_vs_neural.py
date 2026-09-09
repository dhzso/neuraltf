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
    ax.plot(x_grid, kde_a(x_grid), color="#55606B", lw=1.2,
            label=f"All candidates (n = {len(a):,}, med = {np.median(a):.2f})")
    ax.plot(x_grid, kde_n(x_grid), color=C_A, lw=1.5,
            label=f"Neural TFs (n = {len(n):,}, med = {np.median(n):.2f})")

    # Median lines
    ax.axvline(np.median(a), color="#55606B", ls=":", lw=0.9, alpha=0.6)
    ax.axvline(np.median(n), color=C_A, ls=":", lw=0.9, alpha=0.6)

    ks, p = ks_2samp(a, n)
    p_str = "P < 10^{-30}" if p < 1e-30 else f"P = {p:.1e}"

    ax.set_xlabel("Integrated evidence score", fontsize=7.0)
    ax.set_ylabel("Probability density", fontsize=7.0)
    ax.set_title("Integrated Evidence Score Separation: All TFs vs Neural Regulators",
                 fontsize=8.0, fontweight="bold", pad=14)
    ax.text(0.5, 1.02,
            f"Two-sample Kolmogorov–Smirnov test: $D = {ks:.3f}$, ${p_str}$ (Background $n = {len(a):,}$, Neural $n = {len(n):,}$)",
            transform=ax.transAxes, fontsize=6.3, ha="center", va="bottom", color="#444444")
    ax.set_ylim(0, 4.3)

    ax.legend(frameon=False, fontsize=6.2, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save(fig, "03_score_distribution_all_vs_neural")

if __name__=="__main__": build()

