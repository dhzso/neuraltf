"""Bootstrap/Dirichlet confidence intervals on scores for the top-20 candidates.

Reads results/bootstrap_scores_ci.csv (weight-vector uncertainty; see
scripts/stats/bootstrap_confidence.py). Rows are ordered by the OBSERVED
integrated score so the top-20 by real ranking are shown; error bars are
the 95% percentile band from the centered-Dirichlet draw matrix.
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    path = RES / "bootstrap_scores_ci.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing — run scripts/stats/bootstrap_confidence.py first"
        )
    df = pd.read_csv(path)

    # prefer the observed integrated score for ordering (not the bootstrap
    # mean — with weight uncertainty every mean shrinks toward the cohort)
    order_col = "integrated_score" if "integrated_score" in df.columns \
        else ("centered_mean" if "centered_mean" in df.columns else "bootstrap_mean")
    mean_col = "centered_mean" if "centered_mean" in df.columns else "bootstrap_mean"
    lo_col = "centered_ci_95_lo" if "centered_ci_95_lo" in df.columns else "ci_95_lo"
    hi_col = "centered_ci_95_hi" if "centered_ci_95_hi" in df.columns else "ci_95_hi"

    df = df.sort_values(order_col, ascending=True).tail(20)

    fig, ax = plt.subplots(figsize=(8, 6))
    y = np.arange(len(df))
    names = df["gene_name"].values if "gene_name" in df.columns else df["gene_id"].values
    means = df[order_col].values
    lo = df[lo_col].values
    hi = df[hi_col].values

    ax.barh(y, means, height=0.6, color=C_A, alpha=0.7, edgecolor="white", lw=0.3)
    ax.errorbar(means, y, xerr=[np.maximum(means - lo, 0), np.maximum(hi - means, 0)],
                fmt="none", ecolor="#333333", elinewidth=1, capsize=3, capthick=1)

    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8, fontweight="bold")
    ax.set_xlabel("Integrated score (observed; whiskers = 95% band)")
    ax.set_title("Top-20 candidates: 95% uncertainty band under weight resampling\n"
                 "Centered Dirichlet (k=40), 1000 draws — uncertainty in the WEIGHTS, not the genes",
                 fontweight="bold", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "25_bootstrap_ci")

if __name__ == "__main__":
    build()
