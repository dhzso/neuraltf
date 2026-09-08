"""Weight sensitivity — P(Top10) for each candidate."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np

def build():
    sens = load_sens_top10()
    rank_all = load_all()
    top10 = load_top10()
    track_map = dict(zip(top10["gene_id"], top10.get("track",[""]*len(top10))))

    df = sens.sort_values("frac_draws_in_top10", ascending=True)
    fig, ax = plt.subplots(figsize=(W_15COL, 5.2))
    y = np.arange(len(df))
    top10_id_set = set(top10["gene_id"])
    entrant_track = {}
    if "baseline_track" in sens.columns:
        entrant_track = dict(zip(sens["gene_id"], sens["baseline_track"].astype(str)))
    def _color(gid):
        if track_map.get(gid, "") == "A" or entrant_track.get(gid) == "A": return C_A
        if (gid in top10_id_set and track_map.get(gid, "") == "B") or entrant_track.get(gid) == "B": return C_B
        return C_NEURAL
    colors = [_color(r["gene_id"]) for _, r in df.iterrows()]
    bars = ax.barh(y, df["frac_draws_in_top10"].values, color=colors, height=0.62, edgecolor="none")
    
    # Value annotations on top candidates
    for i, (_, r) in enumerate(df.iterrows()):
        val = r["frac_draws_in_top10"]
        if val >= 0.5:
            ax.text(val + 0.015, y[i], f"{val:.1%}", fontsize=6, va="center", color="#333333")

    ax.set_yticks(y)
    ax.set_yticklabels([label(rank_all, g) for g in df["gene_id"]], fontsize=6.5)

    ax.axvline(x=0.8, color="#888888", lw=0.8, ls="--", label="Threshold = 80%")
    ax.set_xlabel("Top 10 retention frequency", fontsize=8)
    ax.set_ylabel("Candidate", fontsize=8)
    ax.set_xlim(0, 1.25)
    ax.set_title("Top 10 retention frequency (1,000 Dirichlet draws)",
                 fontweight="bold", fontsize=8.5, pad=6)
    
    from matplotlib.lines import Line2D
    track_handles = [Line2D([0],[0], marker="s", color="w", markerfacecolor=C_A, markersize=6, label="Track A (benchmark)"),
                     Line2D([0],[0], marker="s", color="w", markerfacecolor=C_B, markersize=6, label="Track B (candidate)"),
                     Line2D([0],[0], marker="s", color="w", markerfacecolor=C_NEURAL, markersize=6, label="Other candidate")]
    handles, labels_leg = ax.get_legend_handles_labels()
    handles.extend(track_handles)
    labels_leg.extend(["Track A (benchmark)","Track B (candidate)","Other candidate"])
    ax.legend(handles, labels_leg, frameon=False, fontsize=6.5, loc="lower right", bbox_to_anchor=(0.98, 0.04))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save(fig, "07_weight_sensitivity_ptop10")

if __name__=="__main__": build()

