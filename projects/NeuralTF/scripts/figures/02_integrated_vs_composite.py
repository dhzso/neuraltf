"""Integrated score vs composite score — showing the effect of domain/GO/ortholog bonuses."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np

def build():
    neural = load_neural()
    top10 = load_top10()
    top10_ids = set(top10["gene_id"].tolist())
    track_map = dict(zip(top10["gene_id"], top10.get("track",[""]*len(top10))))

    fig, ax = plt.subplots(figsize=(6, 5))
    for gid, row in neural.iterrows():
        x = row.get("integrated_score", np.nan)
        y = row.get("composite_score", x)  # fallback if no composite
        if x != y:  # only plot if composite differs from integrated
            color = C_A if track_map.get(row["gene_id"],"")=="A" else C_B if track_map.get(row["gene_id"],"") in ("A","B") else C_NEURAL
            s = 25 if row["gene_id"] in top10_ids else 8
            ax.scatter(x, y, s=s, c=color, alpha=0.6, edgecolors="white", linewidth=0.3, zorder=3)

    # Identity line
    lo, hi = 0, max(neural["integrated_score"].max(), neural.get("composite_score",pd.Series([0])).max())*1.05
    ax.plot([lo,hi],[lo,hi],"--",color="#999999",lw=0.8,label="y = x")

    # Annotate bonus direction
    ax.annotate("Bonuses increase\ncomposite score", xy=(0.55,0.75), fontsize=7,
                color=C_HL, fontweight="bold", ha="center")

    ax.set_xlabel("Integrated evidence score"); ax.set_ylabel("Composite score")
    ax.set_title("Integrated vs composite score\n(bonuses: TF domain +0.05, GO neural +0.03, GO TF +0.02, ortholog +0.02, RNAi brain +0.02)",
                 fontweight="bold", pad=8, fontsize=9)
    ax.legend(frameon=False, fontsize=7); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "02_integrated_vs_composite")

if __name__=="__main__": build()
