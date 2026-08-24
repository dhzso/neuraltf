"""Weight sensitivity — P(Top10) for each candidate."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np

def build():
    sens = load_sens_top10()
    neural = load_neural()
    top10 = load_top10()
    track_map = dict(zip(top10["gene_id"], top10.get("track",[""]*len(top10))))

    df = sens.sort_values("frac_draws_in_top10", ascending=True)
    fig, ax = plt.subplots(figsize=(7, 6))
    y = np.arange(len(df))
    colors = [C_A if track_map.get(r["gene_id"],"")=="A" else C_B if r["gene_id"] in set(top10["gene_id"]) else C_NEURAL for _,r in df.iterrows()]
    bars = ax.barh(y, df["frac_draws_in_top10"].values, color=colors, height=0.65, edgecolor="white", lw=0.3)
    ax.set_yticks(y); ax.set_yticklabels([label(neural, g) for g in df["gene_id"]], fontsize=6)
    for i, (_,r) in enumerate(df.iterrows()):
        c = C_A if track_map.get(r["gene_id"],"")=="A" else C_B if r["gene_id"] in set(top10["gene_id"]) else C_NEURAL
        ax.get_yticklabels()[i].set_color(c)
    ax.axvline(x=0.8, color=C_HL, lw=0.8, ls="--", label="80% threshold")
    ax.set_xlabel("Fraction of 1000 draws in Top 10"); ax.set_xlim(0, 1.05)
    ax.set_title("P(Top 10) under random weight perturbation", fontweight="bold", pad=8)
    ax.legend(frameon=False, fontsize=7); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "07_weight_sensitivity_ptop10")

if __name__=="__main__": build()
