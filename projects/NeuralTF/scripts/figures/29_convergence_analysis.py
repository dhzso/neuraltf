"""Convergence analysis — Dirichlet draw convergence + permutation power curves."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    # Panel A: Convergence of top-10 overlap with increasing draws
    conv_path = RES / "convergence_draws.csv"
    if conv_path.exists():
        conv = pd.read_csv(conv_path)
        draws = conv["n_draws"].values
        overlap = conv["overlap_fraction"].values
        ax1.plot(draws, overlap, color=C_A, lw=2, marker="o", markersize=4)
        ax1.axhline(y=0.9, color="#999999", lw=1, linestyle="--", alpha=0.7, label="90% convergence")
        ax1.set_xlabel("Number of Dirichlet draws")
        ax1.set_ylabel("Fraction of top-10 in final set")
        ax1.set_title("Top-10 candidate overlap converges\nwith increasing draws", fontweight="bold")
        ax1.legend(fontsize=7)
    else:
        ax1.text(0.5, 0.5, "convergence_draws.csv not found", ha="center", va="center",
                 transform=ax1.transAxes, fontsize=9, color="#999999")

    # Panel B: Permutation power curve
    power_path = RES / "permutation_power.csv"
    if power_path.exists():
        power = pd.read_csv(power_path)
        n_perms = power["n_permutations"].values
        power_val = power["power"].values
        ax2.plot(n_perms, power_val, color=C_B, lw=2, marker="s", markersize=4)
        ax2.axhline(y=0.8, color="#999999", lw=1, linestyle="--", alpha=0.7, label="80% power")
        ax2.set_xlabel("Number of permutations")
        ax2.set_ylabel("Statistical power")
        ax2.set_title("Permutation test power increases\nwith more permutations", fontweight="bold")
        ax2.legend(fontsize=7)
    else:
        ax2.text(0.5, 0.5, "permutation_power.csv not found", ha="center", va="center",
                 transform=ax2.transAxes, fontsize=9, color="#999999")

    for ax in [ax1, ax2]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "29_convergence_analysis")

if __name__ == "__main__":
    build()
