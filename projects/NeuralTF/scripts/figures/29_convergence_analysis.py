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
        ax1.set_title("Rank-list convergence with increasing draws\n(real centered-Dirichlet draw matrix, rank-vector Spearman)",
                      fontweight="bold")
        ax1.legend(fontsize=7)
        ax1.set_ylim(0.5, 1.02)
    else:
        ax1.text(0.5, 0.5, "convergence_draws.csv not found\n(run stats/power_analysis.py)",
                 ha="center", va="center", transform=ax1.transAxes,
                 fontsize=9, color="#999999")

    # Panel B: permutation resolution (add-one p-value granularity)
    if res_path.exists():
        res = pd.read_csv(res_path)
        n_perm = res["n_perm"].values
        minp = res["min_detectable_p"].values
        ax2.plot(n_perm, minp, color=C_B, lw=2, marker="s", markersize=4,
                 label="min detectable p = 1/(n+1)")
        ax2.axhline(y=0.05, color=C_HL, lw=1, linestyle="--", alpha=0.9,
                    label="alpha = 0.05")
        ax2.set_xlabel("Number of permutations")
        ax2.set_ylabel("Minimum detectable p-value")
        ax2.set_title("Permutation test p-value granularity vs n_perm\n(add-one estimator; n>=20 needed to reject at 0.05)",
                      fontweight="bold")
        ax2.legend(fontsize=7)
        ax2.set_yscale("log")
    else:
        ax2.text(0.5, 0.5, "permutation_resolution.csv not found\n(run stats/power_analysis.py)",
                 ha="center", va="center", transform=ax2.transAxes,
                 fontsize=9, color="#999999")

    for ax in [ax1, ax2]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout()
    save(fig, "29_convergence_analysis")

if __name__ == "__main__":
    build()
