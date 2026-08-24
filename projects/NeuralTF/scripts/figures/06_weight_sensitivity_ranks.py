"""Weight sensitivity — rank distributions from 1000 random weight draws."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd

def build():
    draws = load_sens_draws()
    sens = load_sens_top10()
    neural = load_neural()
    top10 = load_top10()
    top10_ids = set(top10["gene_id"].tolist())
    track_map = dict(zip(top10["gene_id"], top10.get("track",[""]*len(top10))))

    # Get baseline ranks for all candidates in draws
    candidates = draws["gene_id"].unique()
    baseline = {}
    for gid in candidates:
        row = neural[neural["gene_id"]==gid]
        if len(row)>0:
            baseline[gid] = row.iloc[0].get("integrated_score", 0)
    baseline_ranks = pd.Series(baseline).rank(ascending=False).to_dict()

    fig, ax = plt.subplots(figsize=(8, 5.5))
    y_labels = []
    y_pos = []
    colors = []
    for i, gid in enumerate(sorted(candidates, key=lambda g: baseline_ranks.get(g, 999))):
        sub = draws[draws["gene_id"]==gid]
        if sub.empty: continue
        ranks = sub["rank"].values
        color = C_A if track_map.get(gid,"")=="A" else C_B if gid in top10_ids else C_NEURAL
        alpha = 0.7 if gid in top10_ids else 0.25
        bp = ax.boxplot(ranks, vert=False, positions=[i], widths=0.6,
                        patch_artist=True, showfliers=False,
                        boxprops=dict(facecolor=color, alpha=alpha, edgecolor="none"),
                        medianprops=dict(color="#333", lw=1.2),
                        whiskerprops=dict(color="#999", lw=0.5),
                        capprops=dict(color="#999", lw=0.5))
        y_labels.append(label(neural, gid))
        y_pos.append(i)
        colors.append(color)

    ax.set_yticks(y_pos); ax.set_yticklabels(y_labels, fontsize=5.5)
    for i, c in enumerate(colors):
        ax.get_yticklabels()[i].set_color(c)
    ax.axvline(x=10, color=C_HL, lw=0.8, ls="--", label="Top 10 cutoff")
    ax.set_xlabel("Rank across 1000 weight draws"); ax.invert_yaxis()
    ax.set_title("Weight sensitivity — rank distributions", fontweight="bold", pad=8)
    ax.legend(frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "06_weight_sensitivity_ranks")

if __name__=="__main__": build()
