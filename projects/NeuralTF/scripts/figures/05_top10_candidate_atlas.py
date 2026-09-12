"""Top 10 prioritized candidates — comprehensive publication-grade Candidate Atlas.

Dual-panel 500 DPI figure integrating:
- Panel a: Multi-stream evidence profile matrix (6 streams with numerical values & heatmap)
- Panel b: Prioritization scores (composite bar vs base marker) paired with InterPro TF families & human orthologs
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

TF_FAMILIES = {
    "Homeobox": "#0072B2",
    "bHLH": "#E69F00",
    "Zinc finger": "#009E73",
    "Forkhead": "#D55E00",
    "T-box": "#CC79A7",
    "p53": "#56B4E9",
}


def tf_family_color(domains_str):
    d = str(domains_str).lower()
    if "homeo" in d:
        return "Homeobox", TF_FAMILIES["Homeobox"]
    if "bhlh" in d or "ash" in d:
        return "bHLH", TF_FAMILIES["bHLH"]
    if "znf" in d or "zinc" in d or "c2h2" in d:
        return "Zinc finger", TF_FAMILIES["Zinc finger"]
    if "fork_head" in d or "forkhead" in d or "fox" in d:
        return "Forkhead", TF_FAMILIES["Forkhead"]
    if "t-box" in d or "tbox" in d:
        return "T-box", TF_FAMILIES["T-box"]
    if "p53" in d:
        return "p53", TF_FAMILIES["p53"]
    return "Other", "#888888"


def clean_ortholog_symbol(orth_str):
    s = str(orth_str).strip()
    if not s or s.lower() in ("nan", "none", "-", ""):
        return "—"
    mapping = {
        "forkhead box protein l1": "FOXL1",
        "homeobox protein aristaless-like 4": "ALX4",
        "homeobox protein six6": "SIX6",
        "homeobox protein unc-4 homolog": "UNC4",
        "achaete-scute homolog 3": "ASCL3",
        "pancreas transcription factor 1 subunit alpha": "PTF1A",
        "nk1 transcription factor-related protein 1": "NKX1-1",
        "barh-like 2 homeobox protein": "BARHL2",
        "homeobox protein dlx-1": "DLX1",
        "zinc finger e-box-binding homeobox 1": "ZEB1",
    }
    s_clean = s.split("[")[0].strip().lower()
    for k, v in mapping.items():
        if k in s_clean:
            return v
    words = s.split()[0:3]
    return " ".join(words)


def build():
    neural = load_neural()
    top10 = load_top10()

    if top10.empty:
        fig, ax = plt.subplots(figsize=(7.5, 4.0))
        ax.text(0.5, 0.5, "No candidates", ha="center", va="center")
        save(fig, "05_top10_candidate_atlas")
        return

    # Extract complete multi-dimensional metadata
    records = []
    for _, row in top10.iterrows():
        gid = row["gene_id"]
        track = row.get("track", "B")
        rank = row.get("rank", 1)
        nm = label(neural, gid)
        n_row = neural[neural["gene_id"] == gid]
        base = n_row.iloc[0].get("integrated_score", np.nan) if len(n_row) > 0 else np.nan
        comp = row.get("composite_score", np.nan)
        orth_raw = str(row.get("human_ortholog", "") or "")
        orth_sym = clean_ortholog_symbol(orth_raw)
        domains = str(row.get("interpro_domains", row.get("domains_all", "")))
        fam, fam_c = tf_family_color(domains)

        # Stream values
        s_expr = n_row.iloc[0].get("expression", np.nan) if len(n_row) > 0 else np.nan
        s_spec = n_row.iloc[0].get("specificity", np.nan) if len(n_row) > 0 else np.nan
        s_nspec = n_row.iloc[0].get("neural_specificity", np.nan) if len(n_row) > 0 else np.nan
        s_repro = n_row.iloc[0].get("reproducibility", np.nan) if len(n_row) > 0 else np.nan
        s_lineage = n_row.iloc[0].get("perez_lineage", np.nan) if len(n_row) > 0 else np.nan
        s_inf = n_row.iloc[0].get("perez_influence", np.nan) if len(n_row) > 0 else np.nan

        records.append({
            "gene_id": gid,
            "name": nm,
            "track": track,
            "rank": rank,
            "base": base,
            "composite": comp,
            "ortholog": orth_sym,
            "family": fam,
            "family_color": fam_c,
            "s_expr": s_expr,
            "s_spec": s_spec,
            "s_nspec": s_nspec,
            "s_repro": s_repro,
            "s_lineage": s_lineage,
            "s_inf": s_inf,
        })
    df = pd.DataFrame(records)

    # Separate Track A and Track B, ordered 1..5 top to bottom
    sub_a = df[df["track"] == "A"].sort_values("rank", ascending=True)
    sub_b = df[df["track"] == "B"].sort_values("rank", ascending=True)
    df_sorted = pd.concat([sub_a, sub_b], ignore_index=True)

    n_total = len(df_sorted)
    streams = ["s_expr", "s_spec", "s_nspec", "s_repro", "s_lineage", "s_inf"]
    stream_names = [
        "Expression",
        "Specificity",
        "Neural\nspecificity",
        "Reproduc-\nibility",
        "Perez\nlineage",
        "Perez\ninfluence",
    ]
    mat = df_sorted[streams].values.astype(float)
    masked_mat = np.ma.masked_invalid(mat)

    # Colormap for evidence matrix
    cmap_colors = ["#EDF4F9", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C", "#08306B"]
    base_cmap = mcolors.LinearSegmentedColormap.from_list("pub_blues", cmap_colors)
    base_cmap.set_bad(color="#ECEFF1")

    # Figure setup: 8.6 x 4.8 inches at 500 DPI
    fig = plt.figure(figsize=(8.6, 4.8), dpi=500)

    # Deterministic pixel-perfect axes layout
    # [left, bottom, width, height]
    ax_track = fig.add_axes([0.045, 0.15, 0.015, 0.68])
    ax_mat   = fig.add_axes([0.145, 0.15, 0.355, 0.68])
    ax_score = fig.add_axes([0.550, 0.15, 0.410, 0.68], sharey=ax_mat)

    div_y = 4.5  # Track A is 0..4, Track B is 5..9
    colors = [C_A] * len(sub_a) + [C_B] * len(sub_b)

    # ------------------ TRACK STRIP (Far Left) ------------------
    for i, c in enumerate(colors):
        ax_track.add_patch(plt.Rectangle((0, i - 0.5), 1, 1, color=c, ec="none"))
    ax_track.set_xlim(0, 1)
    ax_track.set_ylim(n_total - 0.5, -0.5)
    ax_track.axhline(div_y, color="white", lw=2.5)
    ax_track.axis("off")

    # Track labels to the left of the track strip ("RNAi+" = screened in
    # King 2024 mmc5; phenotype status is NOT implied — see groundtruth.py)
    fig.text(0.025, 0.65, "Tested (RNAi-screened)\u2020", rotation=90, va="center", ha="center",
             fontsize=6.8, fontweight="bold", color=C_A)
    fig.text(0.025, 0.32, "Not tested", rotation=90, va="center", ha="center",
             fontsize=6.8, fontweight="bold", color=C_B)

    # ------------------ PANEL A: Evidence Stream Matrix ------------------
    panel_tag(ax_mat, "a", x=-0.18, y=1.12)
    im = ax_mat.imshow(masked_mat, aspect="auto", cmap=base_cmap, vmin=0, vmax=1, interpolation="nearest")
    ax_mat.set_xticks(range(len(streams)))
    ax_mat.set_xticklabels(stream_names, fontsize=6.5, fontweight="bold", va="bottom")
    ax_mat.xaxis.tick_top()
    ax_mat.xaxis.set_label_position("top")
    ax_mat.set_yticks(range(n_total))
    ax_mat.set_yticklabels(df_sorted["name"], fontsize=7.2, fontweight="bold")
    ax_mat.set_ylim(n_total - 0.5, -0.5)
    ax_mat.axhline(div_y, color="white", lw=2.5)

    # Display numeric values inside each matrix cell
    for r in range(n_total):
        for c in range(len(streams)):
            val = mat[r, c]
            if np.isnan(val):
                ax_mat.text(c, r, "—", ha="center", va="center", fontsize=5.8, color="#888888")
            else:
                txt_col = "white" if val >= 0.70 else "#222222"
                ax_mat.text(c, r, f"{val:.2f}", ha="center", va="center", fontsize=5.8, fontweight="bold", color=txt_col)

    for s in ax_mat.spines.values():
        s.set_visible(False)
    ax_mat.tick_params(axis="both", length=0, pad=6)

    # ------------------ PANEL B: Prioritization Scores & Annotations ------------------
    panel_tag(ax_score, "b", x=-0.06, y=1.12)
    y_pos = np.arange(n_total)

    # Two-tone stacked horizontal bar chart: Base score + Composite bonus
    base_scores = df_sorted["base"].values
    comp_scores = df_sorted["composite"].values
    bonus_scores = comp_scores - base_scores
    C_BONUS = "#D9822B"  # Distinct warm amber gold for composite bonus

    ax_score.barh(y_pos, base_scores, height=0.56, color=colors, alpha=0.90, edgecolor="none")
    ax_score.barh(y_pos, bonus_scores, left=base_scores, height=0.56, color=C_BONUS, alpha=0.90, edgecolor="none")

    # Column headers above Panel B
    x_score_col = 1.15
    x_fam_col = 1.38
    x_orth_col = 1.62

    ax_score.text(x_score_col, -0.75, "Score", fontsize=6.5, fontweight="bold", ha="center", va="bottom", color="#444444")
    ax_score.text(x_fam_col, -0.75, "TF Family", fontsize=6.5, fontweight="bold", ha="center", va="bottom", color="#444444")
    ax_score.text(x_orth_col, -0.75, "Human Ortholog", fontsize=6.5, fontweight="bold", ha="left", va="bottom", color="#444444")

    # Annotations aligned in neat vertical columns
    for y, (_, r) in enumerate(df_sorted.iterrows()):
        # 1. Composite score
        ax_score.text(x_score_col, y, f"{r['composite']:.2f}", va="center", ha="center", fontsize=6.2, fontweight="bold", color="#222222")
        # 2. TF family badge
        ax_score.text(x_fam_col, y, r['family'], va="center", ha="center", fontsize=6.2,
                      color=r['family_color'], fontweight="bold")
        # 3. Human ortholog (left-aligned)
        ax_score.text(x_orth_col, y, f"Hs: {r['ortholog']}", va="center", ha="left", fontsize=6.2,
                       fontstyle="italic", color="#333333")

    ax_score.axhline(div_y, color="#CCCCCC", lw=1.0, ls="--")
    ax_score.set_xlim(0, 2.05)
    ax_score.set_ylim(n_total - 0.5, -0.5)
    ax_score.set_xlabel("Prioritization Score", fontsize=7.2, fontweight="bold")
    ax_score.set_xticks([0.0, 0.5, 1.0])
    ax_score.tick_params(axis="x", labelsize=6.8)
    ax_score.tick_params(left=False, labelleft=False)
    ax_score.spines["top"].set_visible(False)
    ax_score.spines["right"].set_visible(False)
    ax_score.spines["left"].set_color("#DDDDDD")

    # Colorbar below Panel A
    cbar_ax = fig.add_axes([0.145, 0.038, 0.18, 0.016])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Stream Score", fontsize=6.2, fontweight="bold")
    cbar.set_ticks([0.0, 0.5, 1.0])
    cbar.ax.tick_params(labelsize=6)
    cbar.outline.set_linewidth(0.5)

    # Missing / N/A swatch below Panel A
    fig.patches.append(plt.Rectangle((0.345, 0.038), 0.012, 0.016, transform=fig.transFigure,
                                     facecolor="#ECEFF1", edgecolor="#CCCCCC", lw=0.5, clip_on=False))
    fig.text(0.365, 0.044, "N/A", fontsize=6.0, va="center", color="#555555")

    # Unified Legend below Panel B
    leg_handles = [
        Patch(facecolor=C_A, label="Tested: base score"),
        Patch(facecolor=C_B, label="Not tested: base score"),
        Patch(facecolor=C_BONUS, label="Composite bonus"),
    ]
    ax_score.legend(
        handles=leg_handles,
        loc="lower left",
        bbox_to_anchor=(0.0, -0.22),
        ncol=3,
        frameon=False,
        fontsize=6.2,
    )

    fig.suptitle(
        "Candidate Atlas: High-Confidence & Novel Neural Transcription Factor Regulators",
        fontsize=8.5,
        fontweight="bold",
        y=0.98,
    )

    save(fig, "05_top10_candidate_atlas")


if __name__ == "__main__":
    build()

