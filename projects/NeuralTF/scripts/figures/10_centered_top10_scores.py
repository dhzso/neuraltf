"""Centered Dirichlet — evidence heatmap for top 10 candidates."""
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
    for i, (_, row) in enumerate(order.iterrows()):
        gid = row["gene_id"]
        track = row.get("track","")
        tc = C_A if track=="A" else C_B
        sym = "o" if track=="A" else "s"
        composite = row.get("composite_score", np.nan)
        dirichlet = row.get("dirichlet_median_score", np.nan)
        nm = label(neural, gid)
        ax.barh(y[i]-0.15, composite, height=0.3, color=tc, alpha=0.8, label="Composite" if i==0 else "")
        ax.barh(y[i]+0.15, dirichlet, height=0.3, color=tc, alpha=0.4, label="Dirichlet median" if i==0 else "")
        ax.text(max(composite, dirichlet)+0.005, y[i], nm, fontsize=7, va="center", color=tc, fontweight="bold")
    ax.set_yticks(y); ax.set_yticklabels([f"A" if r["track"]=="A" else "B" for _,r in order.iterrows()], fontsize=7)
    ax.set_xlabel("Score"); ax.set_xlim(0, 1.05)
    ax.set_title("Centered Dirichlet — composite vs median Dirichlet score (Top 10)", fontweight="bold", pad=8)
    ax.legend(frameon=False, fontsize=7, loc="lower right"); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "10_centered_top10_scores")

if __name__=="__main__": build()
