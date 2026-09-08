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
    neural = load_neural()
    proof_map = dict(zip(neural["gene_id"], neural.get("proof_status", [""] * len(neural))))

    fig, ax = plt.subplots(figsize=(W_15COL, 4.4))
    y = np.arange(len(df))
    names = [label(neural, gid) for gid in df["gene_id"]]
    means = df[order_col].values
    lo = df[lo_col].values
    hi = df[hi_col].values

    # Determine track color for each candidate
    colors = []
    for gid in df["gene_id"]:
        ps = str(proof_map.get(gid, "")).lower()
        colors.append(C_A if "validated" in ps or "fstf" in ps else C_B)

    # Point-range plot
    for i in range(len(df)):
        c = colors[i]
        # 95% uncertainty interval
        ax.plot([lo[i], hi[i]], [y[i], y[i]], color=c, lw=1.6, alpha=0.75, zorder=3)
        # Point estimate
        ax.scatter([means[i]], [y[i]], color=c, s=36, edgecolors="white", lw=0.6, zorder=4)

    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=7.5)

    ax.set_xlabel("Integrated score (point = median, bar = 95% CI)", fontsize=8)
    ax.set_ylabel("Candidate", fontsize=8)
    ax.set_title("Score credible intervals (1,000 Dirichlet draws)",
                 fontweight="bold", fontsize=8.5, pad=6)
    
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], color=C_A, marker="o", lw=1.5, markersize=5, label="Track A (benchmark)"),
        Line2D([0], [0], color=C_B, marker="o", lw=1.5, markersize=5, label="Track B (candidate)")
    ]
    ax.legend(handles=legend_handles, loc="upper left", frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "25_bootstrap_ci")

if __name__ == "__main__":
    build()

