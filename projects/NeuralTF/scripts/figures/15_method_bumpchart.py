"""3-method comparison — grouped bar chart of ranks for top-10 candidates."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd

def build():
    fixed_top10 = load_top10()
    centered_top10 = load_centered()
    uniform_top10 = load_uniform()
    neural = load_neural()

    # Build rank lookup from each top-10 file (within-track ranks 1-5)
    def rank_lookup(df):
        d = {}
        for _, row in df.iterrows():
            gid = row["gene_id"]
            d[gid] = row.get("rank", np.nan)
        return d

    fixed_r = rank_lookup(fixed_top10)
    centered_r = rank_lookup(centered_top10)
    uniform_r = rank_lookup(uniform_top10)

    # Use fixed top-10 as the candidate list
    records = []
    for _, row in fixed_top10.iterrows():
        gid = row["gene_id"]
        nm = label(neural, gid)
        track = row.get("track","")
        records.append({
            "name": nm, "track": track,
            "fixed": fixed_r.get(gid, np.nan),
            "centered": centered_r.get(gid, np.nan),
            "uniform": uniform_r.get(gid, np.nan),
        })
    df = pd.DataFrame(records)
    df = df.sort_values("fixed", ascending=True)
    y = np.arange(len(df))
    bw = 0.25

    fig, ax = plt.subplots(figsize=(9, 6))
    colors = {"fixed": C_FIXED, "centered": C_CENTERED, "uniform": C_UNIFORM}
    offsets = {"fixed": -bw, "centered": 0, "uniform": bw}
    methods = ["fixed", "centered", "uniform"]
    labels = ["Fixed-weight", "Centered Dirichlet", "Uniform Dirichlet"]

    for method, lbl in zip(methods, labels):
        vals = df[method].values
        bars = ax.barh(y + offsets[method], vals, height=bw*0.9, color=colors[method],
                       alpha=0.85, edgecolor="white", lw=0.3, label=lbl)
        for i, v in enumerate(vals):
            if pd.notna(v):
                ax.text(v + 0.1, y[i] + offsets[method], f'{int(v)}',
                        fontsize=6, va="center", color=colors[method], fontweight="bold")

    # Gene names with track color
    ax.set_yticks(y)
    ax.set_yticklabels(df["name"], fontsize=8, fontweight="bold")
    for i, track in enumerate(df["track"]):
        ax.get_yticklabels()[i].set_color(C_A if track=="A" else C_B)

    from matplotlib.lines import Line2D
    track_handles = [Line2D([0],[0], marker="s", color="w", markerfacecolor=C_A, markersize=8, label="Track A (RNAi)"),
                     Line2D([0],[0], marker="s", color="w", markerfacecolor=C_B, markersize=8, label="Track B (novel)")]
    handles, labels_leg = ax.get_legend_handles_labels()
    handles.extend(track_handles)
    labels_leg.extend(["Track A (RNAi)", "Track B (novel)"])
    ax.legend(handles, labels_leg, loc="lower right", fontsize=7, frameon=True)

    ax.set_xlabel("Rank within track (1 = highest)", fontsize=9)
    ax.set_ylabel("Candidate", fontsize=9)
    ax.set_title("Rank comparison across three weighting methods (Top 10)\n"
                 "Track A ranks 1-5, Track B ranks 1-5; gene names colored by track",
                 fontweight="bold", pad=10, fontsize=10)
    ax.set_xlim(0, 7)
    ax.set_xticks([1, 2, 3, 4, 5, 6])
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    save(fig, "15_method_bumpchart")

if __name__=="__main__": build()
