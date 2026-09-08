"""Convergence analysis — Dirichlet draw convergence + permutation resolution.

Inputs (written by scripts/stats/power_analysis.py from the REAL draw
matrix and observed scores):
  results/convergence_draws.csv     — n_draws, spearman_vs_full, spearman_std
  results/permutation_resolution.csv — n_perm, min_detectable_p, resolves_p05

Panel B shows the permutation test's p-value GRANULARITY (1/(n+1)) rather
than the retired tautological "power" panel (which was 1.0 by
construction under the label-permutation null).
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    conv_path = RES / "convergence_draws.csv"
    res_path = RES / "permutation_resolution.csv"
    if not conv_path.exists() and not res_path.exists():
        raise FileNotFoundError(
            "Neither convergence_draws.csv nor permutation_resolution.csv found - "
            "run scripts/stats/power_analysis.py (it computes both from the "
            "real Dirichlet draw matrix)."
        )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_2COL, 3.2))

    # Panel A: rank stability vs number of draws (real draw matrix)
    if conv_path.exists():
        conv = pd.read_csv(conv_path)
        x = conv["n_draws"].values
        y = conv["spearman_vs_full"].values
        err = conv["spearman_std"].values if "spearman_std" in conv.columns else None
        if err is not None:
            ax1.errorbar(x, y, yerr=err, color=C_A, lw=1.5, marker="o",
                         markersize=4, capsize=2.5, ecolor="#555555",
                         label="Spearman $r_s$ vs full draws")
        else:
            ax1.plot(x, y, color=C_A, lw=1.5, marker="o", markersize=4,
                     label="Spearman $r_s$ vs full draws")
        ax1.axhline(y=0.95, color="#888888", lw=0.8, linestyle="--", alpha=0.8,
                    label="0.95 stability criterion")
        ax1.set_xlabel("Number of Dirichlet draws", fontsize=8)
        ax1.set_ylabel("Spearman rank stability ($r_s$)", fontsize=8)
        ax1.set_title("Dirichlet draw convergence (k = 40)",
                      fontweight="bold", fontsize=8.5, pad=6)
        ax1.legend(fontsize=6.8, frameon=False, loc="lower right")
        ax1.set_ylim(0.5, 1.03)
    else:
        ax1.text(0.5, 0.5, "convergence_draws.csv not found",
                 ha="center", va="center", transform=ax1.transAxes,
                 fontsize=8, color="#999999")

    # Panel B: permutation resolution (add-one p-value granularity)
    if res_path.exists():
        res = pd.read_csv(res_path)
        n_perm = res["n_perm"].values
        minp = res["min_detectable_p"].values
        ax2.plot(n_perm, minp, color=C_B, lw=1.5, marker="s", markersize=4,
                 label="Granularity: $p_{min} = 1/(n+1)$")
        ax2.axhline(y=0.05, color=C_HL, lw=0.8, linestyle="--", alpha=0.9,
                    label=r"$\alpha = 0.05$ threshold")
        ax2.set_xlabel("Number of permutations ($n$)", fontsize=8)
        ax2.set_ylabel("Minimum detectable $p$-value", fontsize=8)
        ax2.set_title("Permutation test $p$-value resolution floor",
                      fontweight="bold", fontsize=8.5, pad=6)
        ax2.legend(fontsize=6.8, frameon=False, loc="upper right")
        ax2.set_yscale("log")
    else:
        ax2.text(0.5, 0.5, "permutation_resolution.csv not found",
                 ha="center", va="center", transform=ax2.transAxes,
                 fontsize=8, color="#999999")

    for ax in [ax1, ax2]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "29_convergence_analysis")


if __name__ == "__main__":
    build()
