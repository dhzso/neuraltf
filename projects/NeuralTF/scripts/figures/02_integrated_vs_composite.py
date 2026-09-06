"""Integrated score distribution — all 249 TF candidates."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np

def build():
    all249 = load_all()
    scores = all249["integrated_score"].dropna().values

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(scores, bins=30, color=C_A, edgecolor="white", lw=0.4, alpha=0.75)

    med = np.median(scores)
    q25, q75 = np.percentile(scores, [25, 75])
    ax.axvline(med, color=C_B, ls="--", lw=1.2, label=f"Median = {med:.3f}")
    ax.axvline(q25, color="#999", ls=":", lw=0.8, label=f"Q1 = {q25:.3f}")
    ax.axvline(q75, color="#999", ls=":", lw=0.8, label=f"Q3 = {q75:.3f}")

    ax.set_xlabel("Integrated evidence score")
    ax.set_ylabel("Number of TF candidates")
    ax.set_title(f"Integrated scores distribution across TF candidates (n={len(scores)}, median = {med:.3f})",
                 fontweight="bold", pad=8)

    ax.legend(frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "02_integrated_vs_composite")

if __name__=="__main__": build()
