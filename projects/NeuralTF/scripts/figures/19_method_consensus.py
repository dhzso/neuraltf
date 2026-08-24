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

    # Classify each candidate
    all_ids = f_ids | c_ids | u_ids
    records = []
    for gid in all_ids:
        nm = label(neural, gid)
        in_f = gid in f_ids
        in_c = gid in c_ids
        in_u = gid in u_ids
        n_methods = sum([in_f, in_c, in_u])
        track = "A" if gid in f_ids and fixed[fixed["gene_id"]==gid].iloc[0].get("track","") == "A" else "B"
        methods = []
        if in_f: methods.append("Fixed")
        if in_c: methods.append("Centered")
        if in_u: methods.append("Uniform")
        records.append({
            "name": nm, "gene_id": gid, "track": track,
            "n_methods": n_methods,
            "fixed": in_f, "centered": in_c, "uniform": in_u,
            "methods": "+".join(methods),
        })
    df = pd.DataFrame(records).sort_values(["n_methods", "name"], ascending=[False, True])

    # Count by consensus level
    consensus_counts = df["n_methods"].value_counts().sort_index()

    # Consensus color map
    cons_colors = {1: C_B, 2: "#CC79A7", 3: C_A}
    cons_labels = {1: "1 method (method-specific)", 2: "2 methods (partial consensus)", 3: "All 3 methods (robust)"}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1, 1.5]})

    # Left panel: bar chart of consensus counts
    ax = axes[0]
    levels = sorted(consensus_counts.index)
    bars = ax.bar(levels, [consensus_counts[l] for l in levels],
                  color=[cons_colors[l] for l in levels], edgecolor="white", lw=0.5)
    for bar, l in zip(bars, levels):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
                f'{consensus_counts[l]}', ha="center", fontsize=9, fontweight="bold")
    ax.set_xlabel("Number of methods in top-10")
    ax.set_ylabel("Number of candidates")
    ax.set_title("Method consensus", fontweight="bold", pad=8)
    ax.set_xticks(levels)
    ax.set_xticklabels([f"{l}/3" for l in levels])
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    legend_handles = [Patch(facecolor=cons_colors[l], label=cons_labels[l]) for l in levels]
    ax.legend(handles=legend_handles, loc="upper left", fontsize=6, frameon=True)

    # Right panel: table showing which methods each candidate appears in
    ax2 = axes[1]
    ax2.axis("off")
    y = np.arange(len(df))
    for i, (_, r) in enumerate(df.iterrows()):
        yi = len(df) - 1 - i
        # Gene name
        tc = C_A if r["track"] == "A" else C_B
        ax2.text(0.0, yi, r["name"], fontsize=7, va="center", fontweight="bold", color=tc,
                 transform=ax2.transAxes)
        # Method indicators
        for j, (col, lbl) in enumerate([("fixed","F"), ("centered","C"), ("uniform","U")]):
            if r[col]:
                ax2.plot(0.25 + j*0.12, yi, "o", color=cons_colors[r["n_methods"]],
                         markersize=8, transform=ax2.transAxes, zorder=3)
                ax2.text(0.25 + j*0.12, yi, lbl, fontsize=5, ha="center", va="center",
                         color="white", fontweight="bold", transform=ax2.transAxes, zorder=4)
            else:
                ax2.plot(0.25 + j*0.12, yi, "o", color="#DDD", markersize=8,
                         transform=ax2.transAxes, zorder=3)

    # Column headers
    for j, lbl in enumerate(["Fixed", "Centered", "Uniform"]):
        ax2.text(0.25 + j*0.12, len(df) + 0.3, lbl, fontsize=7, ha="center",
                 fontweight="bold", transform=ax2.transAxes)
    ax2.text(0.0, len(df) + 0.3, "Candidate", fontsize=7, fontweight="bold", transform=ax2.transAxes)
    ax2.set_title("Candidate presence across methods", fontweight="bold", pad=8, fontsize=10)

    fig.suptitle("Method consensus analysis — ranking robustness", fontweight="bold", fontsize=11, y=1.02)
    fig.tight_layout(w_pad=2)
    save(fig, "19_method_consensus")

if __name__=="__main__": build()
