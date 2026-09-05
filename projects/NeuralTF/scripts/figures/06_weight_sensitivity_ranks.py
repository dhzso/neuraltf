"""Weight sensitivity — rank distributions from 1000 random weight draws."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd

def build():
    draws = load_sens_draws()
    sens = load_sens_top10()
    rank_all = load_all()
    top10 = load_top10()
    top10_ids = set(top10["gene_id"].tolist())
    track_map = dict(zip(top10["gene_id"], top10.get("track",[""]*len(top10))))

    # Only candidates that ever reach a top-10 slot are drawn (sens already
    # filters to entrants); draws CSV holds shortlist/top-30 rows only.
    sens_ids = set(sens["gene_id"])
    draws = draws[draws["gene_id"].isin(sens_ids)]

    # Baseline ranks from the FULL universe (rank.csv), not the neural view
    candidates = draws["gene_id"].unique()
    base = rank_all.set_index("gene_id")["integrated_score"]
    baseline_ranks = base.rank(ascending=False).to_dict()

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
        y_labels.append(label(rank_all, gid))
        y_pos.append(i)
        colors.append(color)

    ax.set_yticks(y_pos); ax.set_yticklabels(y_labels, fontsize=5.5)
    for i, c in enumerate(colors):
        ax.get_yticklabels()[i].set_color(c)
    ax.axvline(x=10, color=C_HL, lw=0.8, ls="--", label="Top 10 cutoff")
    ax.set_xlabel("Rank across 1000 weight draws (full candidate universe)")
    ax.set_ylabel("TF candidate (sorted by baseline rank)")
    ax.invert_yaxis()
    ax.set_title("Most top-10 candidates maintain stable ranks under weight perturbation",
                 fontweight="bold", pad=8)
    from matplotlib.lines import Line2D
    track_handles = [Line2D([0],[0], marker="s", color="w", markerfacecolor=C_A, markersize=7, label="Track A"),
                     Line2D([0],[0], marker="s", color="w", markerfacecolor=C_B, markersize=7, label="Track B"),
                     Line2D([0],[0], marker="s", color="w", markerfacecolor=C_NEURAL, markersize=7, label="Other entrant")]
    handles, labels_leg = ax.get_legend_handles_labels()
    handles.extend(track_handles)
    labels_leg.extend(["Track A","Track B","Other entrant"])
    ax.legend(handles, labels_leg, frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "06_weight_sensitivity_ranks")

if __name__=="__main__": build()
