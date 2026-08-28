"""Permutation null distribution with real scores overlaid."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    df = pd.read_csv(RES / "permutation_pvalues_full.csv")
    null_scores = df["null_score"].dropna().values
    real_scores = df["real_score"].dropna().values if "real_score" in df.columns else None

    fig, ax = plt.subplots(figsize=(7, 4.5))

    ax.hist(null_scores, bins=50, density=True, color=C_NEURAL, alpha=0.6, edgecolor="white",
            label="Null distribution (permutation)")

    if real_scores is not None:
        for rs in real_scores:
            ax.axvline(x=rs, color=C_HL, lw=1.5, linestyle="--", alpha=0.8)
        ax.axvline(x=real_scores[0], color=C_HL, lw=1.5, linestyle="--",
                   label="Real scores (observed)")

    p_empirical = df["p_empirical"].iloc[0] if "p_empirical" in df.columns else None
    if p_empirical is not None:
        ax.text(0.95, 0.95, f"Empirical p = {p_empirical:.4f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=9,
                fontweight="bold", color=C_HL,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=C_HL, alpha=0.8))

    ax.set_xlabel("Score")
    ax.set_ylabel("Density")
    ax.set_title("Permutation test: real scores significantly exceed null distribution",
                 fontweight="bold", pad=8)
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "26_permutation_null")

if __name__ == "__main__":
    build()
