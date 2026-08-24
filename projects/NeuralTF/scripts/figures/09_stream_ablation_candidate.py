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
    pivot = pivot.reindex(columns=[s for s in STREAM_COLS if s in pivot.columns])
    pivot["_track"] = pivot.index.map(lambda g: track_map.get(g,""))
    pivot = pivot.sort_values(["_track","gene_id"], ascending=[True,True]).drop(columns=["_track"])

    vmax = max(abs(pivot.values.min()), abs(pivot.values.max()), 1)
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(7, 4))
    im = ax.imshow(pivot.values, aspect="auto", cmap=plt.cm.RdBu_r, norm=norm, interpolation="nearest")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([STREAM_L[s] for s in pivot.columns], rotation=45, ha="right", fontsize=8)
    ax.set_xlabel("Evidence stream removed")
    ax.xaxis.set_label_position("top")
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels([label(neural, g) for g in pivot.index], fontsize=7)
    ax.set_ylabel("Candidate")
    for i, gid in enumerate(pivot.index):
        c = C_A if track_map.get(gid,"")=="A" else C_B
        ax.get_yticklabels()[i].set_color(c)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i,j]
            if not np.isnan(v):
                tc = "white" if abs(v)>vmax*0.6 else "#333"
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6, color=tc)
    ax.set_title("Candidate sensitivity — rank change per stream removed", fontweight="bold", pad=12)
    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label("Rank change (− = improves)", fontsize=8)
    from matplotlib.lines import Line2D
    track_handles = [Line2D([0],[0], marker="o", color="w", markerfacecolor=C_A, markersize=6, label="Track A (RNAi)"),
                     Line2D([0],[0], marker="o", color="w", markerfacecolor=C_B, markersize=6, label="Track B (novel)")]
    ax.legend(handles=track_handles, loc="upper right", fontsize=6, frameon=True, title="Gene label color", title_fontsize=7)
    fig.tight_layout(); save(fig, "09_stream_ablation_candidate")

if __name__=="__main__": build()
