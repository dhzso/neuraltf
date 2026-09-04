"""Convergence analysis — Dirichlet draw convergence + permutation power curves.

Inputs (written by scripts/stats/power_analysis.py from the REAL draw
matrix and observed scores):
  results/convergence_draws.csv   — n_draws, spearman_vs_full, spearman_std
  results/permutation_power.csv   — n_perm, power
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    conv_path = RES / "convergence_draws.csv"
    power_path = RES / "permutation_power.csv"
    if not conv_path.exists() and not power_path.exists():
        raise FileNotFoundError(
            "Neither convergence_draws.csv nor permutation_power.csv found — "
            "run scripts/stats/power_analysis.py (it computes both from the "
            "real Dirichlet draw matrix)."
        )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    # Panel A: rank stability vs number of draws (real draw matrix)
    if conv_path.exists():
        conv = pd.read_csv(conv_path)
        x = conv["n_draws"].values
        y = conv["spearman_vs_full"].values
        err = conv["spearman_std"].values if "spearman_std" in conv.columns else None
        if err is not None:
            ax1.errorbar(x, y, yerr=err, color=C_A, lw=2, marker="o",
                         markersize=4, capsize=3, label="Spearman vs full draws")
        else:
            ax1.plot(x, y, color=C_A, lw=2, marker="o", markersize=4,
                     label="Spearman vs full draws")
        ax1.axhline(y=0.95, color="#999999", lw=1, linestyle="--", alpha=0.7,
                    label="0.95 stability")
        ax1.set_xlabel("Number of Dirichlet draws")
        ax1.set_ylabel("Spearman rank stability")
        ax1.set_title("Rank-list convergence with increasing draws\n(real centered-Dirichlet draw matrix)",
                      fontweight="bold")
        ax1.legend(fontsize=7)
        ax1.set_ylim(0.5, 1.02)
    else:
        ax1.text(0.5, 0.5, "convergence_draws.csv not found\n(run stats/power_analysis.py)",
                 ha="center", va="center", transform=ax1.transAxes,
                 fontsize=9, color="#999999")

    # Panel B: permutation power curve (real score vector)
    if power_path.exists():
        power = pd.read_csv(power_path)
        n_perm_col = "n_perm" if "n_perm" in power.columns else "n_permutations"
        ax2.plot(power[n_perm_col].values, power["power"].values,
                 color=C_B, lw=2, marker="s", markersize=4)
        ax2.axhline(y=0.8, color="#999999", lw=1, linestyle="--", alpha=0.7,
                    label="80% power")
        ax2.set_xlabel("Number of permutations")
        ax2.set_ylabel("Estimated power")
        ax2.set_title("Permutation test power vs n_perm\n(observed score distribution, add-one p estimator)",
                      fontweight="bold")
        ax2.legend(fontsize=7)
        ax2.set_ylim(0, 1.05)
    else:
        ax2.text(0.5, 0.5, "permutation_power.csv not found\n(run stats/power_analysis.py)",
                 ha="center", va="center", transform=ax2.transAxes,
                 fontsize=9, color="#999999")

    for ax in [ax1, ax2]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "29_convergence_analysis")

if __name__ == "__main__":
    build()
