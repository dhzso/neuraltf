"""Evidence heatmap for all neural-filtered candidates × 9 streams."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, matplotlib.colors as mcolors

def build():
    neural = load_neural()
    col = "integrated_score"
    
    # Sort with Track A (validated) first, then Track B (novel), then other
    def track_sort_key(r):
        ps = str(r.get("proof_status", "")).lower()
        is_val = 0 if ("validated" in ps or "fstf" in ps) else 1
        return (is_val, -r.get("integrated_score", 0))

    sorted_indices = sorted(neural.index, key=lambda idx: track_sort_key(neural.loc[idx]))
    df = neural.loc[sorted_indices].copy().reset_index(drop=True)
    
    streams = [s for s in STREAM_COLS if s in df.columns]
    mat = df[streams].fillna(0).values

    # Focus on the top 45 neural candidates for maximum readability, or full cohort
    n_show = min(40, len(df))
    df_sub = df.iloc[:n_show]
    mat_sub = mat[:n_show]

    fig = plt.figure(figsize=(W_15COL, 6.2))
    # Layout: [track bar (0.04), heatmap (0.82), colorbar (0.04)]
    gs = fig.add_gridspec(1, 3, width_ratios=[0.05, 0.90, 0.05], wspace=0.08)
    ax_track = fig.add_subplot(gs[0, 0])
    ax_main = fig.add_subplot(gs[0, 1])
    ax_cbar = fig.add_subplot(gs[0, 2])

    # 1. Track sidebar
    track_colors = []
    for _, r in df_sub.iterrows():
        ps = str(r.get("proof_status", "")).lower()
        if "validated" in ps or "fstf" in ps:
            track_colors.append(C_A)
        else:
            track_colors.append(C_B)
    
    ax_track.imshow([[1] for _ in range(n_show)], aspect="auto", cmap="binary", vmin=0, vmax=1)
    for i, color in enumerate(track_colors):
        ax_track.add_patch(plt.Rectangle((-0.5, i - 0.5), 1, 1, color=color, ec="none"))
    ax_track.set_xticks([])
    ax_track.set_yticks([])
    ax_track.set_xlim(-0.5, 0.5)
    ax_track.set_ylim(n_show - 0.5, -0.5)
    ax_track.set_ylabel("Track", fontsize=7.5, fontweight="bold")
    ax_track.spines[:].set_visible(False)

    # 2. Main heatmap
    cmap = plt.cm.YlGnBu
    cmap.set_bad("#F9F9F9")
    im = ax_main.imshow(mat_sub, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    ax_main.set_xticks(range(len(streams)))
    ax_main.set_xticklabels([STREAM_L[s] for s in streams], rotation=40, ha="left", fontsize=7)
    ax_main.xaxis.tick_top()
    ax_main.xaxis.set_label_position("top")
    
    ylabels = [f"{label(neural, gid)} ({df_sub.iloc[i]['integrated_score']:.2f})" 
               for i, gid in enumerate(df_sub["gene_id"])]
    ax_main.set_yticks(range(n_show))
    ax_main.set_yticklabels(ylabels, fontsize=6)
    for i, c in enumerate(track_colors):
        ax_main.get_yticklabels()[i].set_color(c)
        if c == C_A:
            ax_main.get_yticklabels()[i].set_fontweight("bold")
            
    ax_main.set_ylabel("Top neural candidates (sorted by evidence)", fontsize=8)
    ax_main.spines[:].set_visible(False)

    # 3. Colorbar
    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("Evidence score (0–1)", fontsize=7.5)
    cbar.ax.tick_params(labelsize=6.5)

    # Legend for track
    from matplotlib.patches import Patch
    leg_handles = [
        Patch(facecolor=C_A, label="Track A (RNAi-validated)"),
        Patch(facecolor=C_B, label="Track B (novel candidate)")
    ]
    ax_main.legend(handles=leg_handles, loc="upper right", bbox_to_anchor=(1.0, -0.02),
                  ncol=2, frameon=False, fontsize=7)

    fig.suptitle("Evidence stream profiles across prioritized neural TF candidates",
                 fontweight="bold", fontsize=8.5, y=0.99)
    save(fig, "04_evidence_heatmap_neural")

if __name__=="__main__": build()

