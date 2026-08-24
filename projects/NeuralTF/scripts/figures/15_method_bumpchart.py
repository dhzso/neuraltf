"""3-method comparison — slope chart showing rank changes across methods."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd

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
        for name, df in [("fixed",fixed),("centered",centered),("uniform",uniform)]:
            row = df[df["gene_id"]==gid]
            if len(row)>0:
                rec[f"{name}_rank"] = row.iloc[0].get("rank", np.nan)
                rec[f"{name}_score"] = row.iloc[0].get("composite_score", np.nan)
        records.append(rec)
    rank_df = pd.DataFrame(records)
    rank_df = rank_df.dropna(subset=["fixed_rank"], how="all")

    fig, axes = plt.subplots(1, 3, figsize=(10, 5), sharey=False)
    comparisons = [
        ("fixed_rank", "centered_rank", "Fixed", "Centered"),
        ("fixed_rank", "uniform_rank", "Fixed", "Uniform"),
        ("centered_rank", "uniform_rank", "Centered", "Uniform"),
    ]
    for ax, (col_l, col_r, lbl_l, lbl_r) in zip(axes, comparisons):
        sub = rank_df.dropna(subset=[col_l, col_r])
        for _, row in sub.iterrows():
            gid = row["gene_id"]
            track = track_map.get(gid,"")
            color = C_A if track=="A" else C_B
            alpha = 0.9 if gid in set(fixed["gene_id"]) else 0.4
            lw = 1.5 if gid in set(fixed["gene_id"]) else 0.8
            rl, rr = row[col_l], row[col_r]
            ax.plot([0, 1], [rl, rr], "-o", color=color, alpha=alpha, lw=lw, markersize=4, zorder=3)
            if gid in set(fixed["gene_id"]):
                nm = label(neural, gid)
                ax.text(1.05, rr, nm, fontsize=6, va="center", color=color, fontweight="bold")
        ax.set_xlim(-0.1, 1.6)
        ax.set_xticks([0, 1]); ax.set_xticklabels([lbl_l, lbl_r], fontsize=8)
        ax.set_ylabel("Rank"); ax.invert_yaxis()
        ax.set_ylim(10.5, 0.5)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.set_title(f"{lbl_l} vs {lbl_r}", fontsize=9, fontweight="bold")

    fig.suptitle("Rank changes across three weighting methods", fontweight="bold", fontsize=10, y=1.02)
    fig.tight_layout(w_pad=1.5); save(fig, "15_method_bumpchart")

if __name__=="__main__": build()
