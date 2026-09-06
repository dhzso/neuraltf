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

    df = pd.DataFrame(rows).sort_values("base", ascending=True)
    y = np.arange(len(df))

    bonus_cols = ["GO neural", "GO TF", "Human ortholog"]
    bonus_colors = {"GO neural": "#E69F00", "GO TF": "#009E73",
                    "Human ortholog": "#CC79A7"}

    fig, ax = plt.subplots(figsize=(9, 6))

    # Base score bars
    ax.barh(y, df["base"], height=0.6, color="#CCCCCC", edgecolor="white", lw=0.3,
            label="Base integrated score")

    # Stacked bonus bars
    left = df["base"].values.copy()
    for bc in bonus_cols:
        vals = df[bc].values
        mask = vals > 0
        if mask.any():
            ax.barh(y[mask], vals[mask], height=0.6, left=left[mask],
                    color=bonus_colors[bc], alpha=0.85, edgecolor="white", lw=0.3,
                    label=bc)
            left = left + vals

    # Composite score label at end
    for i, (_, r) in enumerate(df.iterrows()):
        comp = r["base"] + r["total_bonus"]
        ax.text(min(comp + 0.005, 1.14), y[i], f'{comp:.3f}', fontsize=7,
                va="center", fontweight="bold")

    # Gene names
    ax.set_yticks(y)
    ax.set_yticklabels(df["name"], fontsize=8, fontweight="bold")
    for i, track in enumerate(df["track"]):
        ax.get_yticklabels()[i].set_color(C_A if track == "A" else C_B)

    ax.set_xlabel("Composite score (base integrated + annotation bonuses)", fontsize=9)
    ax.set_ylabel("Top-10 TF candidate (sorted by base score)", fontsize=9)
    ax.set_title("Annotation bonuses add up to +0.07 on top of base integrated scores for top-10\n"
                 "Gene names colored by track (blue = Track A RNAi-validated, orange = Track B novel)",
                 fontweight="bold", pad=10, fontsize=10)
    ax.legend(loc="lower right", fontsize=6, frameon=True, title="Score component",
             title_fontsize=7)
    ax.set_xlim(0, max(1.15, float(df["composite"].max()) * 1.05) if "composite" in df.columns else 1.15)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); save(fig, "18_composite_bonus_waterfall")

if __name__ == "__main__":
    build()
