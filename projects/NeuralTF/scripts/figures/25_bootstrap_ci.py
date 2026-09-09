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
    ci_df = pd.read_csv(path)
    top10 = load_top10()
    neural = load_neural()

    # Merge credible intervals for top 10 prioritized candidates
    sub = top10.merge(
        ci_df[["gene_id", "centered_mean", "centered_ci_95_lo", "centered_ci_95_hi"]],
        on="gene_id",
        how="left"
    )
    # Sort ascending so top-ranked candidates appear at the top of the horizontal plot
    sub = sub.iloc[::-1].reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(W_15COL, 3.8))
    y = np.arange(len(sub))
    names = [label(neural, gid) for gid in sub["gene_id"]]
    means = sub["centered_mean"].values
    lo = sub["centered_ci_95_lo"].values
    hi = sub["centered_ci_95_hi"].values
    tracks = sub["track"].values

    colors = [C_A if t == "A" else C_B for t in tracks]

    # Point-range plot
    for i in range(len(sub)):
        c = colors[i]
        ax.plot([lo[i], hi[i]], [y[i], y[i]], color=c, lw=1.5, alpha=0.8, zorder=3)
        ax.scatter([means[i]], [y[i]], color=c, s=28, edgecolors="white", lw=0.6, zorder=4)

    # Track separator (between index 4 and 5 in reversed top10)
    ax.axhline(4.5, color="#888888", lw=0.6, ls="--")

    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=6.8)
    ax.set_xlabel("Dirichlet credible score (point = mean, bar = 95% CI)", fontsize=7.0)
    ax.set_ylabel("Candidate", fontsize=7.0)
    ax.set_title("Score credible intervals across 1,000 Dirichlet draws", fontsize=8.0, pad=6)
    
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], color=C_A, marker="o", lw=1.5, markersize=4.5, label="Track A (benchmark)"),
        Line2D([0], [0], color=C_B, marker="o", lw=1.5, markersize=4.5, label="Track B (candidate)")
    ]
    ax.legend(handles=legend_handles, loc="upper left", frameon=False, fontsize=6.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "25_bootstrap_ci")

if __name__ == "__main__":
    build()

