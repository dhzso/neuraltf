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
    "HMG": "#8A4F8B",
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
    if "hmg" in d:
        return "HMG", TF_FAMILIES["HMG"]
    return "Other", "#888888"


def load_family_catalog():
    """Gene -> annotated TF family fallback from the master TF catalog.

    2026-09-23: the InterPro keyword table above only knows six families, so
    any other family rendered a bare "Other" badge (e.g. dd34144/pan =
    TCF7L2, whose domains are "HMG_box_dom; TCF/LEF"). When the keyword
    match misses, look the v6 gene up in the King mmc4 / Perez MOESM5
    annotations (tf_family_king -> tf_family_perez -> tf_class_perez) and
    show that family name instead.
    """
    path = RES.parent / "data" / "master_tf_catalog.csv"
    if not path.exists():
        return {}
    cat = pd.read_csv(path, usecols=["v6_id", "tf_family_king",
                                     "tf_family_perez", "tf_class_perez"])
    out = {}
    for _, r in cat.iterrows():
        for col in ("tf_family_king", "tf_family_perez", "tf_class_perez"):
            v = r[col]
            if pd.notna(v) and str(v).strip() and str(v).strip().lower() != "nan":
                out[str(r["v6_id"])] = str(v).strip()
                break
    return out


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
    family_catalog = load_family_catalog()

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
        is_pheno = bool(row.get("phenotype_confirmed", False))
        disp_nm = f"{nm}†" if is_pheno else nm
        n_row = neural[neural["gene_id"] == gid]
        base = n_row.iloc[0].get("integrated_score", np.nan) if len(n_row) > 0 else np.nan
        comp = row.get("composite_score", np.nan)
        orth_raw = str(row.get("human_ortholog", "") or "")
        orth_sym = clean_ortholog_symbol(orth_raw)
        domains = str(row.get("interpro_domains", row.get("domains_all", "")))
        fam, fam_c = tf_family_color(domains)
        if fam == "Other":
            # Show the annotated family name instead of a bare "Other".
            fam = family_catalog.get(str(gid), "") or fam
            fam_c = TF_FAMILIES.get(fam, fam_c)

        # Extract all 11 stream values
        stream_vals = {}
        for s in STREAM_COLS:
            v = n_row.iloc[0].get(s, np.nan) if len(n_row) > 0 else np.nan
            stream_vals[s] = v

        rec_dict = {
            "gene_id": gid,
            "name": disp_nm,
            "track": track,
            "rank": rank,
            "base": base,
            "composite": comp,
            "ortholog": orth_sym,
            "family": fam,
            "family_color": fam_c,
            "pheno": is_pheno,
        }
        rec_dict.update(stream_vals)
        records.append(rec_dict)
    df = pd.DataFrame(records)

    # Separate Track A and Track B, ordered 1..5 top to bottom
    sub_a = df[df["track"] == "A"].sort_values("rank", ascending=True)
    sub_b = df[df["track"] == "B"].sort_values("rank", ascending=True)
    df_sorted = pd.concat([sub_a, sub_b], ignore_index=True)

    n_total = len(df_sorted)
    streams = STREAM_COLS
    stream_names = [
        "Expression",
        "Specificity",
        "Reproduc-\nibility",
        "RNAi\nscreen",
        "Corr-\nelation",
        "Neural\nenriched",
        "Neural\nspec.",
        "Perez\nlineage",
        "Perez\ninfluence",
        "Fincher\nbrain",
        "Cui\ntemporal",
    ]
    mat = df_sorted[streams].values.astype(float)
    masked_mat = np.ma.masked_invalid(mat)

    # Colormap for evidence matrix
    cmap_colors = ["#EDF4F9", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C", "#08306B"]
    base_cmap = mcolors.LinearSegmentedColormap.from_list("pub_blues", cmap_colors)
    base_cmap.set_bad(color="#ECEFF1")

    # Figure setup: 9.6 x 4.8 inches at 500 DPI
    fig = plt.figure(figsize=(9.6, 5.4), dpi=500)

    # Deterministic pixel-perfect axes layout accommodating all 11 streams
    ax_track = fig.add_axes([0.035, 0.16, 0.012, 0.63])
    ax_mat   = fig.add_axes([0.115, 0.16, 0.490, 0.63])
    ax_score = fig.add_axes([0.655, 0.16, 0.320, 0.63], sharey=ax_mat)

    div_y = 4.5  # Track A is 0..4, Track B is 5..9
    colors = [C_A] * len(sub_a) + [C_B] * len(sub_b)

    # ------------------ TRACK STRIP (Far Left) ------------------
    for i, c in enumerate(colors):
        ax_track.add_patch(plt.Rectangle((0, i - 0.5), 1, 1, color=c, ec="none"))
    ax_track.set_xlim(0, 1)
    ax_track.set_ylim(n_total - 0.5, -0.5)
    ax_track.axhline(div_y, color="white", lw=2.5)
    ax_track.axis("off")

    # Track labels to the left of the track strip
    fig.text(0.025, 0.63, "Track A: RNAi-screened", rotation=90, va="center", ha="center",
             fontsize=6.8, fontweight="bold", color=C_A)
    fig.text(0.025, 0.32, "Track B: not tested", rotation=90, va="center", ha="center",
             fontsize=6.8, fontweight="bold", color=C_B)

    # ------------------ PANEL A: Evidence Stream Matrix ------------------
    panel_tag(ax_mat, "a", x=-0.14, y=1.06)  # keep clear of the 0.878 legend band
    im = ax_mat.imshow(masked_mat, aspect="auto", cmap=base_cmap, vmin=0, vmax=1, interpolation="nearest")
    ax_mat.set_xticks(range(len(streams)))
    ax_mat.set_xticklabels(stream_names, fontsize=5.4, fontweight="bold", va="bottom")
    ax_mat.xaxis.tick_top()
    ax_mat.xaxis.set_label_position("top")
    ax_mat.set_yticks(range(n_total))
    ax_mat.set_yticklabels(df_sorted["name"], fontsize=6.8, fontweight="bold")
    ax_mat.set_ylim(n_total - 0.5, -0.5)
    ax_mat.axhline(div_y, color="white", lw=2.5)

    # Display numeric values inside each matrix cell
    for r in range(n_total):
        for c in range(len(streams)):
            val = mat[r, c]
            if np.isnan(val):
                ax_mat.text(c, r, "—", ha="center", va="center", fontsize=5.2, color="#888888")
            else:
                txt_col = "white" if val >= 0.70 else "#222222"
                ax_mat.text(c, r, f"{val:.2f}", ha="center", va="center", fontsize=5.0, fontweight="bold", color=txt_col)

    for s in ax_mat.spines.values():
        s.set_visible(False)
    ax_mat.tick_params(axis="both", length=0, pad=5)

    # ------------------ PANEL B: Prioritization Scores & Annotations ------------------
    panel_tag(ax_score, "b", x=-0.06, y=1.06)  # keep clear of the 0.878 legend band
    y_pos = np.arange(n_total)

    # Two-tone stacked horizontal bar chart: Base score + Composite bonus
    base_scores = df_sorted["base"].values
    comp_scores = df_sorted["composite"].values
    bonus_scores = comp_scores - base_scores
    C_BONUS = "#D9822B"  # Distinct warm amber gold for composite bonus

    ax_score.barh(y_pos, base_scores, height=0.56, color=colors, alpha=0.90, edgecolor="none")
    ax_score.barh(y_pos, bonus_scores, left=base_scores, height=0.56, color=C_BONUS, alpha=0.90, edgecolor="none")

    # Column headers above Panel B
    x_score_col = 1.12
    x_fam_col = 1.50
    x_orth_col = 1.80

    ax_score.text(x_score_col, -0.75, "Composite", fontsize=6.5, fontweight="bold", ha="center", va="bottom", color="#444444")
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
    ax_score.set_xlim(0, 2.40)
    ax_score.set_ylim(n_total - 0.5, -0.5)
    ax_score.set_xlabel("Composite score (integrated + bonus)", fontsize=7.2, fontweight="bold")
    ax_score.set_xticks([0.0, 0.5, 1.0])
    ax_score.tick_params(axis="x", labelsize=6.8)
    ax_score.tick_params(left=False, labelleft=False)
    ax_score.spines["top"].set_visible(False)
    ax_score.spines["right"].set_visible(False)
    ax_score.spines["left"].set_color("#DDDDDD")

    # Colorbar below Panel A
    cbar_ax = fig.add_axes([0.115, 0.048, 0.16, 0.016])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Evidence score (0–1)", fontsize=6.2, fontweight="bold")
    cbar.set_ticks([0.0, 0.5, 1.0])
    cbar.ax.tick_params(labelsize=6)
    cbar.outline.set_linewidth(0.5)

    # Missing / N/A swatch below Panel A
    fig.patches.append(plt.Rectangle((0.300, 0.048), 0.012, 0.016, transform=fig.transFigure,
                                     facecolor="#ECEFF1", edgecolor="#CCCCCC", lw=0.5, clip_on=False))
    fig.text(0.318, 0.054, "Unassayed", fontsize=6.0, va="center", color="#555555")

    # Unified Legend below Panel B
    leg_handles = [
        Patch(facecolor=C_A, label="Integrated score — Track A (RNAi-screened)"),
        Patch(facecolor=C_B, label="Integrated score — Track B (not tested)"),
        Patch(facecolor=C_BONUS, label="Composite bonus"),
    ]
    fig.legend(
        handles=leg_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.878),
        ncol=3,
        frameon=False,
        fontsize=6.2,
    )

    # Footnote explaining FISH phenotype confirmation
    fig.text(0.360, 0.008,
             "† FISH-confirmed (King et al., 2024)   ·   — unassayed (not measured)",
             fontsize=5.8, fontstyle="italic", color="#555555")

    title_block(
        fig,
        "Candidate Atlas: Dual-Track Top-10 Neural Transcription Factors",
        "Panel a: 11-stream evidence matrix; panel b: composite score = base + bonus (max +0.07)",
    )

    save(fig, "05_top10_candidate_atlas")


if __name__ == "__main__":
    build()
