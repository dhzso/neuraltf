"""3-method comparison — score density distributions."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np
from scipy.stats import gaussian_kde

def build():
    fixed = load_top10()
    centered = load_centered()
    uniform = load_uniform()
    fig, ax = plt.subplots(figsize=(6, 4.5))
    x_grid = np.linspace(0.4, 1.0, 200)
    for df, c, lbl in [(fixed, C_FIXED, "Fixed-weight"), (centered, C_CENTERED, "Centered Dirichlet"), (uniform, C_UNIFORM, "Uniform Dirichlet")]:
        vals = df["composite_score"].dropna().values
        if len(vals)>1:
            kde = gaussian_kde(vals)
            ax.plot(x_grid, kde(x_grid), color=c, lw=1.5, label=lbl)
            ax.fill_between(x_grid, kde(x_grid), color=c, alpha=0.15)
    ax.set_xlabel("Composite score (base integrated + annotation bonuses)")
    ax.set_ylabel("Probability density of candidates per score bin")
    ax.set_title("Three methods produce similar composite score distributions for their top-10 candidates",
                 fontweight="bold", pad=8)
    ax.legend(frameon=False, fontsize=7); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "16_method_score_density")

if __name__=="__main__": build()
