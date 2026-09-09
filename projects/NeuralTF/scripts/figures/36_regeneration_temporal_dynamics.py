"""Figure 36: Regeneration Temporal Dynamics of Prioritized Neural TFs.

Profiles the longitudinal expression trajectories of prioritized transcription
factors across the planarian head regeneration timecourse (Cui et al. 2023):
0d (intact), 6h, 12h, 24h, 2d, 3d, 5d, 7d post-amputation.

Panel a: Heatmap of standardized temporal expression (z-score per candidate)
         across all 8 regeneration stages, separated by Track A and Track B.
Panel b: Relative fold-induction trajectories (Et / E0d) for early-response
         vs late-differentiation neural regulators.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    cui_path = REPO / "projects" / "NeuralTF" / "data" / "cui_atlas_summary.csv"
    if not cui_path.exists():
        raise FileNotFoundError(f"{cui_path} missing")

    cui = pd.read_csv(cui_path)
    top10 = load_top10()
    neural = load_neural()

    time_cols = [
        "mean_cut0d",
        "mean_cut6h",
        "mean_cut12h",
        "mean_cut24h",
        "mean_cut2d",
        "mean_cut3d",
        "mean_cut5d",
        "mean_cut7d",
    ]
    time_labels = ["0d (intact)", "6h", "12h", "24h", "2d", "3d", "5d", "7d"]

    cui["v6_id"] = cui["gene_id"].apply(
        lambda x: f"dd_Smed_v6_{x}" if not str(x).startswith("dd_") else str(x)
    )
    cui = cui.rename(columns={"gene_id": "cui_gene_id"})
    gid_col = "gene_id" if "gene_id" in top10.columns else "gene_id_v6"
    merged = top10.merge(cui, left_on=gid_col, right_on="v6_id", how="left").drop_duplicates(subset=gid_col)

    # Sort: Track A first (ranks 1-5), then Track B (ranks 1-5)
    merged["rank_num"] = pd.to_numeric(merged["rank"], errors="coerce").fillna(99)
    merged = merged.sort_values(["track", "rank_num"], ascending=[True, True]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(4.0, 3.4))

    # Standardize expression per candidate (z-score across timepoints)
    expr_matrix = merged[time_cols].values.astype(float)
    means = expr_matrix.mean(axis=1, keepdims=True)
    stds = expr_matrix.std(axis=1, keepdims=True)
    stds[stds == 0] = 1.0
    z_matrix = (expr_matrix - means) / stds

    n_genes = len(merged)
    gene_labels = []
    for _, r in merged.iterrows():
        nm = label(neural, r[gid_col])
        trk = r["track"]
        rk = int(r["rank_num"]) if r["rank_num"] < 90 else ""
        gene_labels.append(f"{nm}  [{trk}{rk}]")

    # Diverging colormap from muted steel blue to terracotta
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "div_temporal", ["#2B4C6F", "#7AA2C0", "#F4F6F9", "#D98274", "#B04A3E"], N=256
    )

    im = ax.imshow(
        z_matrix,
        aspect="auto",
        cmap=cmap,
        vmin=-2.0,
        vmax=2.0,
        interpolation="nearest",
    )

    ax.set_xticks(range(len(time_labels)))
    ax.set_xticklabels(time_labels, rotation=25, ha="right", fontsize=6.2)
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(gene_labels, fontsize=6.2)

    # Color y-tick labels by Track membership
    for i, ticklabel in enumerate(ax.get_yticklabels()):
        if merged.iloc[i]["track"] == "A":
            ticklabel.set_color(C_A)
        else:
            ticklabel.set_color(C_B)

    # Divider separating Track A and Track B
    ax.axhline(4.5, color="#222222", lw=1.0, ls="--")

    # Legend above the axes indicating Track membership
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=C_A, edgecolor="#222222", lw=0.5, label="Track A (RNAi validated)"),
        Patch(facecolor=C_B, edgecolor="#222222", lw=0.5, label="Track B (Novel candidates)"),
    ]
    ax.legend(
        handles=legend_elements,
        loc="lower left",
        bbox_to_anchor=(0.0, 1.02),
        ncol=2,
        frameon=False,
        fontsize=6.0,
        handletextpad=0.4,
        columnspacing=1.0,
    )

    ax.set_xlabel("Regeneration stage post-amputation", fontsize=7.0)
    ax.set_title("Regeneration Temporal Dynamics (Cui et al. 2023)", fontsize=8.0, pad=18)

    # Subtle cell borders
    ax.set_xticks(np.arange(-0.5, len(time_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_genes, 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.04)
    cbar.set_label("Standardized expression ($z$-score)", fontsize=6.8)
    cbar.ax.tick_params(labelsize=6.0)

    fig.tight_layout()
    save(fig, "36_regeneration_temporal_dynamics")

if __name__ == "__main__":
    build()
