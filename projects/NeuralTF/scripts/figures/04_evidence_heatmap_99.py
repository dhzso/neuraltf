"""Evidence heatmap for all 99 neural-filtered candidates × 7 streams."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, matplotlib.colors as mcolors

def build():
    neural = load_neural()
    col = "integrated_score"
    if col in neural.columns:
        order = neural.sort_values(col, ascending=False).index
    else:
        order = neural.index
    df = neural.loc[order]
    mat = df[STREAM_COLS].fillna(0).values

    fig, ax = plt.subplots(figsize=(10, 12))
    cmap = plt.cm.YlOrRd; cmap.set_bad("#F0F0F0")
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(len(STREAM_COLS)))
    ax.set_xticklabels([STREAM_L[s] for s in STREAM_COLS], rotation=45, ha="right", fontsize=8)
    ax.set_xlabel("Evidence stream")
    ax.xaxis.set_label_position("top")
    ylabels = [label(neural, gid) for gid in df["gene_id"]]
    ax.set_yticks(range(len(df))); ax.set_yticklabels(ylabels, fontsize=5.5)
    ax.set_ylabel("Candidate")
    for i, gid in enumerate(df["gene_id"]):
        ps = str(neural.iloc[i].get("proof_status",""))
        if "validated" in ps.lower() or "fstf" in ps.lower():
            ax.get_yticklabels()[i].set_color(C_A)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i,j]
            if v > 0:
                tc = "white" if v > 0.6 else "#333"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=4, color=tc)
    ax.set_title("Evidence heatmap — all 99 neural candidates", fontweight="bold", pad=12, fontsize=10)
    cbar = fig.colorbar(im, ax=ax, fraction=0.015, pad=0.02)
    cbar.set_label("Evidence strength", fontsize=8)
    from matplotlib.lines import Line2D
    track_handles = [Line2D([0],[0], marker="o", color="w", markerfacecolor=C_A, markersize=6, label="RNAi-validated"),
                     Line2D([0],[0], marker="o", color="w", markerfacecolor="#333", markersize=6, label="Other")]
    ax.legend(handles=track_handles, loc="upper right", fontsize=6, frameon=True, title="Gene label color", title_fontsize=7)
    fig.tight_layout(); save(fig, "04_evidence_heatmap_99")

if __name__=="__main__": build()
