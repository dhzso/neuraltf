"""Stream ablation — candidate sensitivity heatmap."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd, matplotlib.colors as mcolors

def _ablate(df, exclude_idx):
    scores = np.zeros(len(df))
    for i, (_, row) in enumerate(df.iterrows()):
        vals = np.array([row.get(s, np.nan) for s in STREAM_COLS], dtype=float)
        w = W.copy(); w[exclude_idx] = 0
        present = ~np.isnan(vals) & (vals != 0)
        if present.any():
            scores[i] = np.sum(w[present]*vals[present]) / np.sum(w[present])
    return scores

def build():
    neural = load_neural()
    top10 = load_top10()
    top10_ids = set(top10["gene_id"].tolist())
    track_map = dict(zip(top10["gene_id"], top10.get("track",[""]*len(top10))))
    base_scores = _ablate(neural, -1)
    base_ranks = pd.Series(base_scores).rank(ascending=False).values

    rows = []
    for i, (_, row) in enumerate(neural.iterrows()):
        gid = row["gene_id"]
        for j, s in enumerate(STREAM_COLS):
            ab = _ablate(neural, j)
            ab_ranks = pd.Series(ab).rank(ascending=False).values
            rows.append({"gene_id":gid, "stream":s, "rank_change":ab_ranks[i]-base_ranks[i]})

    df = pd.DataFrame(rows)
    top10_sub = df[df["gene_id"].isin(top10_ids)]
    pivot = top10_sub.pivot_table(index="gene_id", columns="stream", values="rank_change", aggfunc="first")
    pivot = pivot.reindex(columns=[s for s in STREAM_HEATMAP_ORDER if s in pivot.columns])
    pivot["_track"] = pivot.index.map(lambda g: track_map.get(g,""))
    pivot = pivot.sort_values(["_track","gene_id"], ascending=[True,True]).drop(columns=["_track"])

    vmax = max(abs(pivot.values.min()), abs(pivot.values.max()), 1)
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    fig = plt.figure(figsize=(W_15COL, 3.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.035, 0.925, 0.04], wspace=0.03)
    ax_track = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1])
    ax_cbar = fig.add_subplot(gs[0, 2])

    # 1. Track sidebar with candidate labels on the left
    track_colors = [C_A if i < 5 else C_B for i in range(len(pivot))]
    ax_track.imshow([[1] for _ in range(len(pivot))], aspect="auto", cmap="binary", vmin=0, vmax=1)
    for i, color in enumerate(track_colors):
        ax_track.add_patch(plt.Rectangle((-0.5, i - 0.5), 1, 1, color=color, ec="none"))
    ylabels = [label(neural, g) for g in pivot.index]
    ax_track.set_xticks([])
    ax_track.set_yticks(range(len(pivot)))
    ax_track.set_yticklabels(ylabels, fontsize=6.5)
    ax_track.tick_params(left=False, right=False, length=0)
    ax_track.set_xlim(-0.5, 0.5)
    ax_track.set_ylim(len(pivot) - 0.5, -0.5)
    ax_track.spines[:].set_visible(False)

    # 2. Main heatmap
    im = ax.imshow(pivot.values, aspect="auto", cmap=plt.cm.RdBu_r, norm=norm, interpolation="nearest")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([STREAM_L[s] for s in pivot.columns], rotation=40, ha="left", fontsize=6.5)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.set_yticks([])
    ax.tick_params(length=0)
        
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            if not np.isnan(v):
                tc = "white" if abs(v) > vmax * 0.45 else "#222222"
                sign_str = f"+{int(v)}" if v > 0 else (f"{int(v)}" if v < 0 else "0")
                ax.text(j, i, sign_str, ha="center", va="center", fontsize=5.8, color=tc)

    # Divider line between Track A and Track B
    ax.axhline(4.5, color="#555555", lw=0.8, ls="--")
    ax_track.axhline(4.5, color="#555555", lw=0.8, ls="--")
    ax.spines[:].set_visible(False)
    
    # 3. Colorbar
    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("Rank shift (Δrank)", fontsize=6.8)
    cbar.ax.tick_params(labelsize=6.0)

    # Legend for track
    from matplotlib.patches import Patch
    leg_handles = [
        Patch(facecolor=C_A, label="Track A (benchmark)"),
        Patch(facecolor=C_B, label="Track B (candidate)")
    ]
    ax.legend(handles=leg_handles, loc="upper right", bbox_to_anchor=(1.0, -0.06),
              ncol=2, frameon=False, fontsize=6.2)

    fig.suptitle("Candidate stream sensitivity (Δrank)",
                 fontsize=8.0, y=0.98)
    fig.subplots_adjust(left=0.18, right=0.92, top=0.78, bottom=0.10)
    save(fig, "09_stream_ablation_candidate")

if __name__=="__main__": build()

