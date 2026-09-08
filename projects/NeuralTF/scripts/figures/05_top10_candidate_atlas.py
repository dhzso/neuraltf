"""Top 10 candidates — domain-colored bar chart with track-colored gene names."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# TF family colors (from InterPro domain annotations)
TF_FAMILIES = {
    "Homeobox": "#0072B2",
    "bHLH": "#E69F00",
    "Zinc finger": "#009E73",
    "Fork head": "#D55E00",
    "T-box": "#CC79A7",
    "p53": "#56B4E9",
}

def tf_family_color(domains_str):
    d = str(domains_str)
    for fam, c in TF_FAMILIES.items():
        if fam.lower() in d.lower():
            return fam, c
    return "Other", "#AAAAAA"

def build():
    neural = load_neural()
    top10 = load_top10()

    if top10.empty:
        fig, ax = plt.subplots(figsize=(W_15COL, 4))
        ax.text(0.5, 0.5, "No candidates", ha="center", va="center")
        save(fig, "05_top10_candidate_atlas")
        return

    # Extract clean details
    records = []
    for _, row in top10.iterrows():
        gid = row["gene_id"]
        track = row.get("track", "B")
        nm = label(neural, gid)
        n_row = neural[neural["gene_id"] == gid]
        base = n_row.iloc[0].get("integrated_score", np.nan) if len(n_row) > 0 else np.nan
        comp = row.get("composite_score", np.nan)
        orth = str(row.get("human_ortholog", "") or "")
        orth_clean = orth if orth.lower() not in ("nan", "none", "") else "—"
        domains = str(row.get("interpro_domains", row.get("domains_all", "")))
        fam, _ = tf_family_color(domains)
        records.append({
            "gene_id": gid, "name": nm, "track": track,
            "base": base, "composite": comp,
            "ortholog": orth_clean, "family": fam,
        })
    df = pd.DataFrame(records)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(W_15COL, 4.8), sharex=True,
                                   gridspec_kw={"height_ratios": [1, 1], "hspace": 0.35})

    tracks_data = [("Track A: RNAi-validated benchmark TFs (n = 5)", "A", ax1, C_A),
                   ("Track B: Novel high-confidence targets for validation (n = 5)", "B", ax2, C_B)]

    for title, track_code, ax, color in tracks_data:
        sub = df[df["track"] == track_code].sort_values("composite", ascending=True)
        y = np.arange(len(sub))

        # Horizontal bars for composite score
        bars = ax.barh(y, sub["composite"], height=0.55, color=color, alpha=0.88,
                       edgecolor="none", label="Composite score")
        
        # Overlay point for base score
        ax.scatter(sub["base"], y, color="#212121", s=30, zorder=4,
                   label="Base score (multi-atlas evidence)", edgecolors="white", lw=0.6)

        # Annotations on the right
        for i, (_, r) in enumerate(sub.iterrows()):
            ax.text(r["composite"] + 0.015, y[i], 
                    f"{r['composite']:.3f}  [{r['family']}, {r['ortholog']}]",
                    fontsize=6.5, va="center", color="#333333")

        ax.set_yticks(y)
        ax.set_yticklabels(sub["name"], fontsize=8, fontweight="bold")
        ax.set_title(title, fontweight="bold", fontsize=8, loc="left", pad=5, color=color)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(0, 1.25)

    ax2.set_xlabel("Prioritization score (base evidence + annotation bonuses)", fontsize=8)
    ax1.legend(loc="lower right", frameon=False, fontsize=7)

    fig.suptitle("NeuralTF Top-10 Prioritized Transcription Factors",
                 fontweight="bold", fontsize=9, y=0.98)
    save(fig, "05_top10_candidate_atlas")

if __name__=="__main__": build()

