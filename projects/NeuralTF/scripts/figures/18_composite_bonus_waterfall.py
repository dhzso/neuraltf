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
    s3 = pd.read_csv(RES / "supplementary_table_S3_centered_all_candidates.csv")
    neural = load_neural()
    top10 = load_top10()
    top10_ids = set(top10["gene_id"])
    track_map = dict(zip(top10["gene_id"], top10["track"]))

    sub = s3[s3["gene_id"].isin(top10_ids)].copy()
    rows = []
    for _, row in sub.iterrows():
        gid = row["gene_id"]
        nm = label(neural, gid)
        track = track_map.get(gid, "B")
        base = row["integrated_score"]
        b_neural = row["bonus_go_neural"]
        b_tf = row["bonus_go_tf"]
        b_orth = row["bonus_human_ortholog"]
        comp = base + b_neural + b_tf + b_orth
        rows.append({
            "gene_id": gid,
            "name": nm,
            "track": track,
            "base": base,
            "GO neural": b_neural,
            "GO TF": b_tf,
            "Human ortholog": b_orth,
            "composite": comp
        })

    df = pd.DataFrame(rows)
    # Order by track (A first, then B) and within track by composite descending
    df = df.sort_values(["track", "composite"], ascending=[False, True]).reset_index(drop=True)
    y = np.arange(len(df))

    bonus_cols = ["GO neural", "GO TF", "Human ortholog"]
    bonus_colors = {
        "GO neural": "#B04A3E",       # terracotta
        "GO TF": "#4A7C59",           # sage green
        "Human ortholog": "#5C82A6"   # steel blue
    }

    fig, ax = plt.subplots(figsize=(W_15COL, 3.8))

    # Base score bars
    ax.barh(y, df["base"], height=0.58, color="#788896", edgecolor="none",
            label="Base evidence score")

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

    # Composite score label at end of each bar
    for i, (_, r) in enumerate(df.iterrows()):
        comp = r["composite"]
        ax.text(comp + 0.015, y[i], f"{comp:.3f}", fontsize=6.2,
                va="center", color="#333333")

    # Track dividing line
    ax.axhline(4.5, color="#888888", lw=0.6, ls="--")

    # Gene names
    ax.set_yticks(y)
    ax.set_yticklabels(df["name"], fontsize=6.8)

    ax.set_xlabel("Prioritization score (base score + additive bonuses)", fontsize=7.0)
    ax.set_ylabel("Candidate", fontsize=7.0)
    ax.set_title("Score composition and bonus breakdown (Top 10 TFs)",
                 fontsize=8.0, pad=16)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=4, frameon=False, fontsize=6.2)
    ax.set_xlim(0, 1.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save(fig, "18_composite_bonus_waterfall")


if __name__ == "__main__":
    build()
