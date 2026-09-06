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
    fig, ax = plt.subplots(figsize=(7, 6))
    y = np.arange(len(df))
    top10_id_set = set(top10["gene_id"])
    # 2026-09-06: entrants' actual tracks from the challengers CSV
    entrant_track = {}
    if "baseline_track" in sens.columns:
        entrant_track = dict(zip(sens["gene_id"], sens["baseline_track"].astype(str)))
    def _color(gid):
        if track_map.get(gid, "") == "A" or entrant_track.get(gid) == "A": return C_A
        if (gid in top10_id_set and track_map.get(gid, "") == "B") or entrant_track.get(gid) == "B": return C_B
        return C_NEURAL
    colors = [_color(r["gene_id"]) for _, r in df.iterrows()]
    bars = ax.barh(y, df["frac_draws_in_top10"].values, color=colors, height=0.65, edgecolor="white", lw=0.3)
    ax.set_yticks(y); ax.set_yticklabels([label(rank_all, g) for g in df["gene_id"]], fontsize=6)
    for i, (_,r) in enumerate(df.iterrows()):
        ax.get_yticklabels()[i].set_color(_color(r["gene_id"]))
    # 2026-09-06: the 80% line carries no statistical justification —
    # relabeled as a descriptive reference, and the title no longer
    # claims a threshold result.
    ax.axvline(x=0.8, color=C_HL, lw=0.8, ls="--", label="80% (descriptive reference)")
    ax.set_xlabel("Fraction of 1000 draws in the dual-track Top 10")
    ax.set_ylabel("TF candidate (sorted by P(Top 10))")
    ax.set_xlim(0, 1.05)
    ax.set_title("Candidate robustness: fraction of weight draws retained in the top 10",
                 fontweight="bold", pad=8)
    from matplotlib.lines import Line2D
    track_handles = [Line2D([0],[0], marker="s", color="w", markerfacecolor=C_A, markersize=7, label="Track A"),
                     Line2D([0],[0], marker="s", color="w", markerfacecolor=C_B, markersize=7, label="Track B"),
                     Line2D([0],[0], marker="s", color="w", markerfacecolor=C_NEURAL, markersize=7, label="Other entrant")]
    handles, labels_leg = ax.get_legend_handles_labels()
    handles.extend(track_handles)
    labels_leg.extend(["Track A","Track B","Other entrant"])
    ax.legend(handles, labels_leg, frameon=False, fontsize=7); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "07_weight_sensitivity_ptop10")

if __name__=="__main__": build()
