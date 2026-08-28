"""Bootstrap confidence intervals on integrated scores for top-20 candidates."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    df = pd.read_csv(RES / "bootstrap_scores_ci.csv")
    df = df.sort_values("mean_score", ascending=True).tail(20)

    fig, ax = plt.subplots(figsize=(8, 6))
    y = np.arange(len(df))
    names = df["gene_name"].values if "gene_name" in df.columns else df["gene_id"].values
    means = df["mean_score"].values
    ci_low = df["ci_low"].values
    ci_high = df["ci_high"].values

    errors_low = means - ci_low
    errors_high = ci_high - means

    ax.barh(y, means, height=0.6, color=C_A, alpha=0.7, edgecolor="white", lw=0.3)
    ax.errorbar(means, y, xerr=[errors_low, errors_high], fmt="none", ecolor="#333333",
                elinewidth=1, capsize=3, capthick=1)

    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8, fontweight="bold")
    ax.set_xlabel("Integrated score (mean ± 95% CI)")
    ax.set_title("Top-20 candidates with bootstrap confidence intervals\n"
                 "Error bars = 95% CI from 1000 bootstrap resamples",
                 fontweight="bold", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "25_bootstrap_ci")

if __name__ == "__main__":
    build()
