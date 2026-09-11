"""Weight sensitivity — P(Top10) for each candidate."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np

def build():
    sens = load_sens_top10()
    rank_all = load_all()
    top10 = load_top10()
    track_map = dict(zip(top10["gene_id"], top10.get("track", [""] * len(top10))))

    df = sens.sort_values("frac_draws_in_top10", ascending=True)
    fig, ax = plt.subplots(figsize=(W_15COL, 5.2))
    y = np.arange(len(df))
    top10_id_set = set(top10["gene_id"])
    entrant_track = {}
    if "baseline_track" in sens.columns:
        entrant_track = dict(zip(sens["gene_id"], sens["baseline_track"].astype(str)))

    def _color(gid):
        if track_map.get(gid, "") == "A" or entrant_track.get(gid) == "A":
            return C_A
        if (gid in top10_id_set and track_map.get(gid, "") == "B") or entrant_track.get(gid) == "B":
            return C_B
        return C_NEURAL

    colors = [_color(r["gene_id"]) for _, r in df.iterrows()]
    bars = ax.barh(y, df["frac_draws_in_top10"].values, color=colors, height=0.62, edgecolor="none")

    # Value annotations on all candidates with P >= 20%
    for i, (_, r) in enumerate(df.iterrows()):
        val = r["frac_draws_in_top10"]
        if val >= 0.20:
            ax.text(val + 0.015, y[i], f"{val:.1%}", fontsize=5.8, va="center", color="#333333")

    ax.set_yticks(y)
    ax.set_yticklabels([label(rank_all, g) for g in df["gene_id"]], fontsize=5.8)

    # Reference lines for majority consensus and high confidence
    ax.axvline(x=0.8, color="#888888", lw=0.7, ls=":")
    ax.axvline(x=0.5, color="#555555", lw=0.8, ls="--")

    ax.set_xlabel("Posterior probability of top-10 ranking P(Top 10) under weight perturbations", fontsize=7.0)
    ax.set_ylabel("Candidate and challenger neural TFs (N = 58)", fontsize=7.0)
    ax.set_xlim(0, 1.14)
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels(["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=6.2)

    ax.set_title("Centered Dirichlet Prior Sensitivity: Top-10 Inclusion Probability",
                 fontsize=8.0, fontweight="bold", pad=14)
    ax.text(0.5, 1.02,
            "Empirical retention frequency across 1,000 Dirichlet draws (k=40) for prioritized candidates & challengers (N = 58)",
            transform=ax.transAxes, fontsize=6.2, ha="center", va="bottom", color="#444444")

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [
        Line2D([0], [0], color="#888888", lw=0.7, ls=":", label="High confidence (P ≥ 80%)"),
        Line2D([0], [0], color="#555555", lw=0.8, ls="--", label="Majority consensus (P ≥ 50%)"),
        Patch(facecolor=C_A, edgecolor="none", label="Tested (n = 5)\u2020"),
        Patch(facecolor=C_B, edgecolor="none", label="Not tested (n = 5)"),
        Patch(facecolor=C_NEURAL, edgecolor="none", label="Challenger TF (n = 48)"),
    ]
    ax.legend(handles=handles, frameon=True, facecolor="white", edgecolor="#D0D7DE",
              fontsize=5.8, loc="lower right", bbox_to_anchor=(0.98, 0.03))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save(fig, "07_weight_sensitivity_ptop10")

if __name__=="__main__": build()

