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

    f_ids = set(fixed["gene_id"])
    c_ids = set(centered["gene_id"])
    u_ids = set(uniform["gene_id"])

    all_ids = f_ids | c_ids | u_ids
    records = []
    for gid in all_ids:
        in_f = gid in f_ids
        in_c = gid in c_ids
        in_u = gid in u_ids
        n_methods = sum([in_f, in_c, in_u])
        records.append({"gene_id": gid, "n_methods": n_methods})
    df = pd.DataFrame(records)

    consensus_counts = df["n_methods"].value_counts().sort_index()
    cons_colors = {1: C_B, 2: "#CC79A7", 3: C_A}
    cons_labels = {1: "1 method (method-specific)", 2: "2 methods (partial)", 3: "All 3 methods (robust)"}

    fig, ax = plt.subplots(figsize=(5, 4.5))
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
    ax.legend(handles=legend_handles, loc="upper left", fontsize=7, frameon=True)

    fig.tight_layout()
    save(fig, "19_method_consensus")

if __name__=="__main__": build()
