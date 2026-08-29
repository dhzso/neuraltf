"""Perez influence comparison — TF classification vs influence score distributions."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

def build():
    # Try loading from MOESM files or fallback to results
    moesm5_path = RES / "MOESM5 ES.csv"
    moesm19_path = RES / "MOESM19 ES.csv"

    if moesm5_path.exists() and moesm19_path.exists():
        df_class = pd.read_csv(moesm5_path)
        df_infl = pd.read_csv(moesm19_path)
    else:
        # Fallback: use available data
        df_class = pd.read_csv(RES / "supplementary_table_S1_method_comparison.csv")
        df_infl = pd.read_csv(RUN / "rank.csv")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    has_labels = False

    # Panel A: Classification distributions
    score_col = "composite_score" if "composite_score" in df_class.columns else "integrated_score"
    if "classification" in df_class.columns:
        groups = df_class["classification"].unique()
        for grp in groups:
            subset = df_class[df_class["classification"] == grp]
            if score_col in subset.columns:
                ax1.hist(subset[score_col].dropna(), bins=20, alpha=0.5, label=grp)
                has_labels = True
    elif "track" in df_class.columns:
        for track in df_class["track"].unique():
            subset = df_class[df_class["track"] == track]
            if score_col in subset.columns:
                ax1.hist(subset[score_col].dropna(), bins=15, alpha=0.5, label=f"Track {track}")
                has_labels = True

    ax1.set_xlabel("Score")
    ax1.set_ylabel("Count")
    ax1.set_title("Score distribution by TF classification", fontweight="bold")
    if has_labels:
        ax1.legend(fontsize=7)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: Perez influence score distribution
    if "perez_influence" in df_infl.columns:
        infl = df_infl["perez_influence"].dropna().values
        x_grid = np.linspace(infl.min(), infl.max(), 200)
        if len(infl) > 1:
            kde = gaussian_kde(infl)
            ax2.fill_between(x_grid, kde(x_grid), color=STREAM_C["perez_influence"], alpha=0.3)
            ax2.plot(x_grid, kde(x_grid), color=STREAM_C["perez_influence"], lw=2)
        ax2.axvline(x=np.median(infl), color=C_HL, lw=1.5, linestyle="--",
                    label=f"Median = {np.median(infl):.3f}")
        ax2.set_xlabel("Perez influence score")
        ax2.set_ylabel("Density")
        ax2.set_title("Perez ANANSE influence score distribution", fontweight="bold")
        ax2.legend(fontsize=7)
    else:
        ax2.text(0.5, 0.5, "perez_influence column not found", ha="center", va="center",
                 transform=ax2.transAxes, fontsize=9, color="#999999")

    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "32_perez_influence_comparison")

if __name__ == "__main__":
    build()
