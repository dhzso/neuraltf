"""Method consensus — which candidates appear in top-10 across methods."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from matplotlib.patches import Patch

def build():
    fixed = load_top10()
    centered = load_centered()
    uniform = load_uniform()
    neural = load_neural()

    f_ids = set(fixed["gene_id"])
    c_ids = set(centered["gene_id"])
    u_ids = set(uniform["gene_id"])

    all_ids = f_ids | c_ids | u_ids
    records = []
    for gid in all_ids:
        nm = label(neural, gid)
        in_f = gid in f_ids
        in_c = gid in c_ids
        in_u = gid in u_ids
        n_methods = sum([in_f, in_c, in_u])
        track = "A" if in_f and fixed[fixed["gene_id"]==gid].iloc[0].get("track","") == "A" else "B"
        records.append({
            "name": nm, "gene_id": gid, "track": track,
            "n_methods": n_methods,
            "fixed": in_f, "centered": in_c, "uniform": in_u,
        })
    df = pd.DataFrame(records).sort_values(["n_methods", "name"], ascending=[False, True]).reset_index(drop=True)

    consensus_counts = df["n_methods"].value_counts().sort_index()
    cons_colors = {1: C_B, 2: "#CC79A7", 3: C_A}
    cons_labels = {1: "1 method (method-specific)", 2: "2 methods (partial)", 3: "All 3 methods (robust)"}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), gridspec_kw={"width_ratios": [1, 1.8]})

    # Left panel: bar chart
    ax = axes[0]
    levels = sorted(consensus_counts.index)
    bars = ax.bar(levels, [consensus_counts[l] for l in levels],
                  color=[cons_colors[l] for l in levels], edgecolor="white", lw=0.5, width=0.6)
    for bar, l in zip(bars, levels):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.15,
                f'{consensus_counts[l]}', ha="center", fontsize=10, fontweight="bold")
    ax.set_xlabel("Number of methods in top-10")
    ax.set_ylabel("Number of unique TF candidates")
    ax.set_title("Most candidates appear in only 1 method's top-10,\n"
                 "but the 3-method overlap is substantial", fontweight="bold", pad=8)
    ax.set_xticks(levels)
    ax.set_xticklabels([f"{l}/3" for l in levels])
    ax.set_ylim(0, max(consensus_counts) + 1.5)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    legend_handles = [Patch(facecolor=cons_colors[l], label=cons_labels[l]) for l in levels]
    ax.legend(handles=legend_handles, loc="upper left", fontsize=6, frameon=True)

    # Right panel: presence matrix using imshow
    ax2 = axes[1]
    n = len(df)
    methods_cols = ["fixed", "centered", "uniform"]
    matrix = np.zeros((n, 3))
    for i, (_, r) in enumerate(df.iterrows()):
        for j, col in enumerate(methods_cols):
            matrix[i, j] = 1 if r[col] else 0

    im = ax2.imshow(matrix, cmap="YlOrRd", aspect="auto", vmin=0, vmax=1, interpolation="nearest")
    ax2.set_xticks(range(3))
    ax2.set_xticklabels(["Fixed", "Centered", "Uniform"], fontsize=8, fontweight="bold")
    ax2.xaxis.set_label_position("top")
    ax2.xaxis.tick_top()

    # Y-axis: gene names
    ax2.set_yticks(range(n))
    ax2.set_yticklabels([df.iloc[i]["name"] for i in range(n)], fontsize=7, fontweight="bold")
    for i in range(n):
        tc = C_A if df.iloc[i]["track"] == "A" else C_B
        ax2.get_yticklabels()[i].set_color(tc)

    # Annotate cells
    for i in range(n):
        for j in range(3):
            present = matrix[i, j] == 1
            txt = "●" if present else "○"
            tc = "white" if present else "#CCC"
            ax2.text(j, i, txt, ha="center", va="center", fontsize=12, color=tc)

    # Consensus label on right margin
    for i in range(n):
        nm = df.iloc[i]["n_methods"]
        c = cons_colors[nm]
        ax2.text(3.3, i, f"{nm}/3", fontsize=7, va="center", fontweight="bold", color=c)

    ax2.set_title("Which candidates are supported by 1, 2, or all 3 methods?", fontweight="bold", pad=8, fontsize=10, y=1.05)
    # Legend for right panel
    from matplotlib.lines import Line2D
    dot_handles = [Line2D([0],[0], marker="o", color="w", markerfacecolor="#E69F00", markersize=8, label="Present in method"),
                   Line2D([0],[0], marker="o", color="w", markerfacecolor="#DDD", markersize=8, label="Absent"),
                   Patch(facecolor=C_A, label="3/3 (robust)"),
                   Patch(facecolor="#CC79A7", label="2/3 (partial)"),
                   Patch(facecolor=C_B, label="1/3 (specific)")]
    ax2.legend(handles=dot_handles, loc="lower right", fontsize=5.5, frameon=True, ncol=2)

    fig.suptitle("Three methods agree on 3 candidates, supporting robust prioritization", fontweight="bold", fontsize=11, y=1.02)
    fig.tight_layout(w_pad=2)
    save(fig, "19_method_consensus")

if __name__=="__main__": build()
