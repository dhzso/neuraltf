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

    fig, ax = plt.subplots(figsize=(W_15COL, 3.8))
    im = ax.imshow(pivot.values, aspect="auto", cmap=plt.cm.RdBu_r, norm=norm, interpolation="nearest")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([STREAM_L[s] for s in pivot.columns], rotation=38, ha="left", fontsize=7)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    
    ylabels = [label(neural, g) for g in pivot.index]
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels(ylabels, fontsize=7)
    
    for i, gid in enumerate(pivot.index):
        c = C_A if track_map.get(gid, "") == "A" else C_B
        ax.get_yticklabels()[i].set_color(c)
        ax.get_yticklabels()[i].set_fontweight("bold")
        
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            if not np.isnan(v):
                tc = "white" if abs(v) > vmax * 0.55 else "#222222"
                sign_str = f"+{int(v)}" if v > 0 else (f"{int(v)}" if v < 0 else "0")
                ax.text(j, i, sign_str, ha="center", va="center", fontsize=6.5, color=tc, fontweight="bold")

    # Add divider line between Track A (rows 0-4) and Track B (rows 5-9)
    ax.axhline(4.5, color="#333333", lw=1.2, ls="--")
    ax.text(-0.85, 2.0, "Track A\n(RNAi)", ha="center", va="center", fontsize=6.5, color=C_A, fontweight="bold", rotation=90)
    ax.text(-0.85, 7.0, "Track B\n(Novel)", ha="center", va="center", fontsize=6.5, color=C_B, fontweight="bold", rotation=90)

    ax.set_ylabel("Top-10 candidate", fontsize=8)
    ax.spines[:].set_visible(False)
    
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.03)
    cbar.set_label("Rank shift (+ drops, − improves)", fontsize=7.5)
    cbar.ax.tick_params(labelsize=6.5)

    ax.set_title("Candidate sensitivity to individual stream removal (Δrank upon omission)",
                 fontweight="bold", fontsize=8.5, pad=18)
    fig.tight_layout()
    save(fig, "09_stream_ablation_candidate")

if __name__=="__main__": build()

