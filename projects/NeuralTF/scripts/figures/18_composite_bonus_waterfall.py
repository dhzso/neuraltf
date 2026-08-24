"""Composite bonus waterfall — shows base score + each bonus for top-10."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd

GO_NEURAL_KW = {"neurogenesis","nervous system","neural","brain","axon","dendrite",
                "synapse","neuron","ganglion","central nervous"}
GO_TF_KW = {"transcription","DNA-binding","transcription factor","regulation of gene expression",
            "sequence-specific DNA binding","RNA polymerase II"}

def compute_bonuses(row):
    bonuses = {}
    # TF domain
    dom = str(row.get("dna_binding_domains", row.get("domains_all", "")))
    tf_flag = str(row.get("mmc4_tf_flag", ""))
    bonuses["TF domain"] = 0.05 if (dom.strip() and dom != "nan") or tf_flag == "TF" else 0.0
    # GO neural
    go_terms = str(row.get("go_terms", "")).lower()
    bonuses["GO neural"] = 0.03 if any(kw in go_terms for kw in GO_NEURAL_KW) else 0.0
    # GO TF
    bonuses["GO TF"] = 0.02 if any(kw in go_terms for kw in GO_TF_KW) else 0.0
    # Human ortholog
    orth = str(row.get("human_ortholog", row.get("planmine_human_ortholog_desc", "")))
    bonuses["Human ortholog"] = 0.02 if (orth.strip() and orth != "nan") else 0.0
    # Brain RNAi
    proof = str(row.get("proof_status", ""))
    bonuses["Brain RNAi"] = 0.02 if "validated" in proof.lower() else 0.0
    return bonuses

def build():
    s2 = pd.read_csv(RES / "supplementary_table_S2_fixed_all249.csv")
    neural = load_neural()
    top10 = load_top10()
    top10_ids = set(top10["gene_id"])

    rows = []
    for _, row in s2[s2["gene_id"].isin(top10_ids)].iterrows():
        base = row.get("integrated_score", 0)
        bonuses = compute_bonuses(row)
        total_bonus = sum(bonuses.values())
        nm = label(neural, row["gene_id"])
        track = "A" if row.get("proof_status","") == "known_rnai_validated" else "B"
        rows.append({"name": nm, "track": track, "base": base, **bonuses, "total_bonus": total_bonus})

    df = pd.DataFrame(rows).sort_values("base", ascending=True)
    y = np.arange(len(df))

    bonus_cols = ["TF domain", "GO neural", "GO TF", "Human ortholog", "Brain RNAi"]
    bonus_colors = {"TF domain": "#0072B2", "GO neural": "#E69F00", "GO TF": "#009E73",
                    "Human ortholog": "#CC79A7", "Brain RNAi": "#D55E00"}

    fig, ax = plt.subplots(figsize=(9, 6))

    # Base score bars
    ax.barh(y, df["base"], height=0.6, color="#CCCCCC", edgecolor="white", lw=0.3, label="Base score")

    # Stacked bonus bars
    left = df["base"].values
    for bc in bonus_cols:
        vals = df[bc].values
        mask = vals > 0
        if mask.any():
            ax.barh(y[mask], vals[mask], height=0.6, left=left[mask],
                    color=bonus_colors[bc], alpha=0.85, edgecolor="white", lw=0.3, label=bc)
            left = left + vals

    # Composite score label at end
    for i, (_, r) in enumerate(df.iterrows()):
        comp = r["base"] + r["total_bonus"]
        ax.text(comp + 0.005, y[i], f'{comp:.3f}', fontsize=7, va="center", fontweight="bold")

    # Gene names
    ax.set_yticks(y)
    ax.set_yticklabels(df["name"], fontsize=8, fontweight="bold")
    for i, track in enumerate(df["track"]):
        ax.get_yticklabels()[i].set_color(C_A if track == "A" else C_B)

    ax.set_xlabel("Score", fontsize=9)
    ax.set_ylabel("Candidate", fontsize=9)
    ax.set_title("Composite score breakdown — base score + annotation bonuses (Top 10)\n"
                 "Gene names colored by track (blue = Track A RNAi-validated, orange = Track B novel)",
                 fontweight="bold", pad=10, fontsize=10)
    ax.legend(loc="lower right", fontsize=6, frameon=True, title="Score component", title_fontsize=7)
    ax.set_xlim(0, 1.15)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "18_composite_bonus_waterfall")

if __name__ == "__main__": build()
