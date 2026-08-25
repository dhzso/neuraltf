"""3-method comparison — rank correlation heatmap."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from scipy.stats import spearmanr

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
                rec[f"{name}_score"] = row.iloc[0].get("composite_score", np.nan)
        records.append(rec)
    score_df = pd.DataFrame(records).dropna()

    methods = ["fixed","centered","uniform"]
    labels = ["Fixed","Centered","Uniform"]
    n = len(methods)
    corr = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i==j: corr[i,j]=1.0
            elif i<j:
                rho, _ = spearmanr(score_df[f"{methods[i]}_score"], score_df[f"{methods[j]}_score"])
                corr[i,j] = rho; corr[j,i] = rho

    fig, ax = plt.subplots(figsize=(5, 4.5))
    im = ax.imshow(corr, cmap="RdYlBu_r", vmin=0.8, vmax=1.0, aspect="equal")
    ax.set_xticks(range(n)); ax.set_xticklabels(labels, fontsize=9)
    ax.set_yticks(range(n)); ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Weighting method (column)")
    ax.set_ylabel("Weighting method (row)")
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{corr[i,j]:.3f}", ha="center", va="center", fontsize=10,
                    color="white" if corr[i,j]<0.9 else "#333", fontweight="bold")
    ax.set_title("All three weighting methods produce highly concordant rankings (Spearman rho > 0.9)",
                 fontweight="bold", pad=8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Spearman rho", fontsize=8)
    fig.tight_layout(); save(fig, "17_method_rank_correlation")

if __name__=="__main__": build()
