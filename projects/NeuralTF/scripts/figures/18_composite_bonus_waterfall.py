"""Composite bonus waterfall — shows base score + each bonus for top-10.

Bonus components and weights mirror prioritize.apply_bonuses EXACTLY:
GO neural +0.03, GO TF +0.02, human ortholog +0.02 (total +0.07). The old
version invented five components (TF domain +0.05, Brain RNAi +0.02) that
do not exist in the scoring code.
"""
from __future__ import annotations
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[3] / "src"))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd

from bioforge.projects.neuraltf.planmine import go_term_flags


def compute_bonuses(row):
    bonuses = {}
    go_terms = str(row.get("go_terms", "") or "")
    has_neural_go = has_go_tf = False
    seen = set()
    for term in go_terms.split(";"):
        t = term.strip()
        if not t or t.lower() in ("nan", "none") or t.lower() in seen:
            continue
        seen.add(t.lower())
        is_neural, is_tf = go_term_flags(t)
        has_neural_go = has_neural_go or is_neural
        has_go_tf = has_go_tf or is_tf
    bonuses["GO neural"] = 0.03 if has_neural_go else 0.0
    bonuses["GO TF"] = 0.02 if has_go_tf else 0.0
    orth = str(row.get("human_ortholog",
                       row.get("planmine_human_ortholog_desc", "")) or "")
    bonuses["Human ortholog"] = 0.02 if orth.strip() and orth.lower() != "nan" else 0.0
    return bonuses

def build():
    s2 = pd.read_csv(RES / "supplementary_table_S2_fixed_all_candidates.csv")
    neural = load_neural()
    top10 = load_top10()
    top10_ids = set(top10["gene_id"])

    rows = []
    for _, row in s2[s2["gene_id"].isin(top10_ids)].iterrows():
        base = row.get("integrated_score", 0)
        bonuses = compute_bonuses(row)
        total_bonus = sum(bonuses.values())
        nm = label(neural, row["gene_id"])
        track = "A" if row.get("proof_status", "") == "known_rnai_validated" else "B"
        rows.append({"name": nm, "track": track, "base": base, **bonuses,
                     "total_bonus": total_bonus})

    df = pd.DataFrame(rows)
    df["composite"] = df["base"] + df["total_bonus"]
    df = df.sort_values("composite", ascending=True)
    y = np.arange(len(df))

    bonus_cols = ["GO neural", "GO TF", "Human ortholog"]
    bonus_colors = {"GO neural": "#B04A3E", "GO TF": "#4A7C59",
                    "Human ortholog": "#5C82A6"}

    fig, ax = plt.subplots(figsize=(W_15COL, 4.2))

    # Base score bars
    ax.barh(y, df["base"], height=0.58, color="#788896", edgecolor="none",
            label="Base score")

    # Stacked bonus bars
    left = df["base"].values.copy()
    for bc in bonus_cols:
        vals = df[bc].values
        mask = vals > 0
        if mask.any():
            ax.barh(y[mask], vals[mask], height=0.58, left=left[mask],
                    color=bonus_colors[bc], alpha=0.9, edgecolor="none",
                    label=f"{bc} (+{vals[mask][0]:.2f})")
            left = left + vals

    # Composite score label at end
    for i, (_, r) in enumerate(df.iterrows()):
        comp = r["composite"]
        ax.text(comp + 0.015, y[i], f"{comp:.2f}", fontsize=7,
                va="center", color="#222222")

    # Gene names
    ax.set_yticks(y)
    ax.set_yticklabels(df["name"], fontsize=7.5)

    ax.set_xlabel("Composite score", fontsize=8)
    ax.set_ylabel("Candidate", fontsize=8)
    ax.set_title("Score composition (top 10 TFs)",
                 fontweight="bold", fontsize=8.5, pad=18)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=4, frameon=False, fontsize=6.8)
    ax.set_xlim(0, 1.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save(fig, "18_composite_bonus_waterfall")


if __name__ == "__main__":
    build()
