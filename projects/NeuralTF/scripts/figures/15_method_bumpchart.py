"""3-method comparison — bump chart (rank flow across fixed/centered/uniform)."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from matplotlib.lines import Line2D

def build():
    fixed = load_top10()
    centered = load_centered()
    uniform = load_uniform()
    neural = load_neural()
    track_map = dict(zip(fixed["gene_id"], fixed.get("track",[""]*len(fixed))))

    all_ids = set()
    for df in [fixed, centered, uniform]:
        if "gene_id" in df.columns: all_ids.update(df["gene_id"].tolist())

    records = []
    for gid in all_ids:
        rec = {"gene_id":gid}
        for name, df, sc in [("fixed",fixed,"composite_score"),("centered",centered,"composite_score"),("uniform",uniform,"composite_score")]:
            row = df[df["gene_id"]==gid]
            if len(row)>0:
                rec[f"{name}_rank"] = row.iloc[0].get("rank", np.nan)
        records.append(rec)
    rank_df = pd.DataFrame(records).dropna(subset=["fixed_rank","centered_rank","uniform_rank"], thresh=2)

    methods = ["Fixed", "Centered", "uniform"]
    x = np.arange(len(methods))
    fig, ax = plt.subplots(figsize=(6, 5))
    for _, row in rank_df.iterrows():
        gid = row["gene_id"]
        track = track_map.get(gid,"")
        color = C_A if track=="A" else C_B if track=="B" else C_NEURAL
        ranks = [row.get(f"{m.lower()}_rank", np.nan) for m in methods]
        valid = [(i,r) for i,r in enumerate(ranks) if pd.notna(r)]
        if len(valid)<2: continue
        xv = [v[0] for v in valid]; yv = [v[1] for v in valid]
        ax.plot(xv, yv, "-o", color=color, alpha=0.8, markersize=5, lw=1.2, zorder=3)
        last_x, last_y = xv[-1], yv[-1]
        ax.annotate(label(neural, gid), (last_x, last_y), fontsize=6, xytext=(5,0),
                    textcoords="offset points", va="center", color=color, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(["Fixed\nweight","Dirichlet\ncentered","Dirichlet\nuniform"], fontsize=8)
    ax.set_ylabel("Rank in Top 10"); ax.invert_yaxis(); ax.set_ylim(0.5, 10.5)
    ax.set_title("Rank stability across three weighting methods", fontweight="bold", pad=8)
    legend_elements = [Line2D([0],[0],color=C_A,marker="o",label="Track A",lw=1),
                       Line2D([0],[0],color=C_B,marker="o",label="Track B",lw=1)]
    ax.legend(handles=legend_elements, frameon=False, fontsize=7, loc="lower right")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "15_method_bumpchart")

if __name__=="__main__": build()
