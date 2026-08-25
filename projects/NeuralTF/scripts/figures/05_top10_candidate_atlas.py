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

    records = []
    for _, row in top10.iterrows():
        gid = row["gene_id"]
        track = row.get("track","")
        nm = label(neural, gid)
        n_row = neural[neural["gene_id"]==gid]
        integrated = n_row.iloc[0].get("integrated_score", np.nan) if len(n_row)>0 else np.nan
        comp_fixed = row.get("composite_score", np.nan)
        proof = str(row.get("proof_status",""))
        domains = str(row.get("interpro_domains", row.get("domains_all","")))
        fam, fam_c = tf_family_color(domains)
        records.append({
            "gene_id": gid, "name": nm, "track": track,
            "integrated": integrated, "composite": comp_fixed,
            "proof": proof, "domains": domains,
            "family": fam, "family_color": fam_c,
        })
    df = pd.DataFrame(records)
    df = df.sort_values("composite", ascending=True)
    y = np.arange(len(df))

    fig, ax = plt.subplots(figsize=(9, 6))

    # Bars colored by TF family
    bar_colors = df["family_color"].tolist()
    ax.barh(y, df["composite"], height=0.65, color=bar_colors, alpha=0.8,
            edgecolor="white", lw=0.5)

    # Score annotations
    for i, (_, r) in enumerate(df.iterrows()):
        ax.text(r["composite"] + 0.008, y[i], f'{r["composite"]:.3f}',
                fontsize=7, va="center", fontweight="bold")

    # Gene names (color-coded by track)
    ax.set_yticks(y)
    ax.set_yticklabels(df["name"], fontsize=8, fontweight="bold")
    for i, (_, r) in enumerate(df.iterrows()):
        ax.get_yticklabels()[i].set_color(C_A if r["track"]=="A" else C_B)

    # Proof status annotation
    for i, (_, r) in enumerate(df.iterrows()):
        proof_short = r["proof"].replace("known_rnai_validated","RNAi validated").replace("novel_domain","novel")[:15]
        ax.text(0.005, y[i] + 0.3, proof_short, fontsize=5.5, va="bottom", color="#666", fontstyle="italic")

    # TF family legend
    seen = []
    legend_handles = []
    for _, r in df.iterrows():
        if r["family"] not in [s[0] for s in seen]:
            seen.append((r["family"], r["family_color"]))
    for fam, c in seen:
        if fam != "Other":
            legend_handles.append(Patch(facecolor=c, alpha=0.8, label=fam))
    legend_handles.append(Patch(facecolor="#AAAAAA", alpha=0.8, label="Other TF"))
    ax.legend(handles=legend_handles, loc="lower right", fontsize=6, frameon=True,
              title="TF family (GO/InterPro)", title_fontsize=7)

    # Track legend (for gene name colors)
    track_handles = [Line2D([0],[0], marker="o", color="w", markerfacecolor=C_A, markersize=8, label="Track A (RNAi-validated)"),
                     Line2D([0],[0], marker="o", color="w", markerfacecolor=C_B, markersize=8, label="Track B (novel)")]
    ax2 = ax.twinx()
    ax2.set_yticks([])
    ax2.legend(handles=track_handles, loc="upper right", fontsize=6, frameon=True, title="Track", title_fontsize=7)

    ax.set_xlabel("Composite score (fixed-weight)", fontsize=9)
    ax.set_ylabel("TF candidate (ranked by composite score)", fontsize=9)
    ax.set_title("Top 10 candidates span multiple TF families\n"
                 "Bars colored by TF family; gene names colored by track",
                 fontweight="bold", pad=10, fontsize=10)
    ax.set_xlim(0, max(df["composite"]) + 0.12)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    fig.tight_layout(); save(fig, "05_top10_candidate_atlas")

if __name__=="__main__": build()
