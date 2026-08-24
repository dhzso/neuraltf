"""Centered Dirichlet — method-specific base score vs composite for top 10."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np

def build():
    centered = load_centered()
    neural = load_neural()
    order = centered.sort_values("composite_score", ascending=False)
    fig, ax = plt.subplots(figsize=(7, 5))
    y = np.arange(len(order))[::-1]
    names = []
    for i, (_, row) in enumerate(order.iterrows()):
        gid = row["gene_id"]
        track = row.get("track","")
        tc = C_A if track=="A" else C_B
        base = row.get("dirichlet_median_score", np.nan)
        comp = row.get("composite_score", np.nan)
        nm = label(neural, gid)
        names.append(nm)
        ax.barh(y[i]-0.15, comp, height=0.3, color=tc, alpha=0.85, edgecolor="white", lw=0.3)
        ax.barh(y[i]+0.15, base, height=0.3, color=tc, alpha=0.35, edgecolor="white", lw=0.3)
        ax.text(max(comp, base)+0.005, y[i], nm, fontsize=7, va="center", color=tc, fontweight="bold")
        # Show bonus
        bonus = comp - base if pd.notna(comp) and pd.notna(base) else 0
        ax.text(comp + 0.02, y[i] - 0.15, f"+{bonus:.3f}", fontsize=5, va="center", color="#666")
    ax.set_yticks(y); ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel("Score"); ax.set_xlim(0, 1.05)
    ax.set_title("Centered Dirichlet (k=40) — Dirichlet median vs composite (Top 10)\n"
                 "Dark = composite (median + bonuses), Light = Dirichlet median only",
                 fontweight="bold", pad=8, fontsize=9)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "10_centered_top10_scores")

if __name__=="__main__": build()
