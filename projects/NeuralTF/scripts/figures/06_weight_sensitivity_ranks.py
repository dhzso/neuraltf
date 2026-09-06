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

    # 2026-09-06 audit fix: the draws CSV now persists FULL 1000-draw
    # rank vectors for every entrant/top-30 gene (the old per-draw
    # persistence filter kept only the good draws for oscillating genes,
    # biasing these boxplots toward stability — e.g. Zeb-1's median read
    # 49 vs the true 76). An in-script completeness check enforces it.
    sens_ids = set(sens["gene_id"])
    draws = draws[draws["gene_id"].isin(sens_ids)]
    n_draws = draws["draw"].max()
    counts = draws.groupby("gene_id").size()
    bad = counts[counts != n_draws]
    if len(bad) > 0:
        raise AssertionError(
            f"fig 06: truncated rank vectors for {len(bad)} genes "
            f"(e.g. {dict(list(bad.items())[:3])}) — regenerate "
            f"weight_sensitivity_draws.csv with the fixed "
            f"run_weight_sensitivity.py")

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
        # 2026-09-06: entrants' actual tracks come from the challengers
        # CSV (baseline_track column), not the 10-gene baseline map —
        # RNAi-validated entrants like dd12722 were mislabeled "Other".
        entrant_track = ""
        if "baseline_track" in sens.columns:
            m = sens.loc[sens["gene_id"]==gid, "baseline_track"]
            if len(m): entrant_track = str(m.iloc[0])
        color = C_A if track_map.get(gid,"")== "A" or entrant_track=="A" \
            else C_B if (gid in top10_ids and track_map.get(gid,"")=="B") or entrant_track=="B" \
            else C_NEURAL
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
    # 2026-09-06: the reference line is the DUAL-TRACK SHORTLIST boundary
    # (top 5 Track A + top 5 Track B), not full-universe rank 10 — the
    # two "10"s are different quantities and must not be conflated.
    ax.axvline(x=30, color=C_HL, lw=0.8, ls="--",
               label="Top-30 full-universe rank (shortlist entrant zone)")
    ax.set_xlabel("Rank across 1000 weight draws (full candidate universe; complete per-gene rank vectors)")
    ax.set_ylabel("TF candidate (sorted by baseline rank)")
    ax.invert_yaxis()
    ax.set_title("Rank stability under weight perturbation (complete draw distributions)",
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
