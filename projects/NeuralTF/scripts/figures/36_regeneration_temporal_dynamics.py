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

    time_cols = ["mean_cut0d", "mean_cut6h", "mean_cut12h", "mean_cut24h",
                 "mean_cut2d", "mean_cut3d", "mean_cut5d", "mean_cut7d"]
    time_labels = ["0d (intact)", "6h", "12h", "24h", "2d", "3d", "5d", "7d"]

    cui["v6_id"] = cui["gene_id"].apply(
        lambda x: f"dd_Smed_v6_{x}" if not str(x).startswith("dd_") else str(x)
    )
    cui = cui.rename(columns={"gene_id": "cui_raw_id"})
    gid_col = "gene_id" if "gene_id" in top10.columns else "gene_id_v6"
    merged = top10.merge(cui, left_on=gid_col, right_on="v6_id", how="left").drop_duplicates(subset=gid_col)

    # Sort: Track A first, then Track B, each by composite score descending
    score_col = "composite_score" if "composite_score" in merged.columns else "integrated_score"
    merged = merged.sort_values(["track", score_col], ascending=[True, False]).reset_index(drop=True)

    fig = plt.figure(figsize=(W_2COL, 3.4))
    gs = fig.add_gridspec(1, 4, width_ratios=[0.03, 0.43, 0.025, 0.49], wspace=0.16)
    ax_track = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[0, 1])
    ax_cbar = fig.add_subplot(gs[0, 2])
    ax_traj = fig.add_subplot(gs[0, 3])

    # --- Panel a: Temporal Expression Heatmap ---
    expr_matrix = merged[time_cols].values.astype(float)
    # Standardize per gene (z-score) to visualize temporal kinetics
    means = expr_matrix.mean(axis=1, keepdims=True)
    stds = expr_matrix.std(axis=1, keepdims=True)
    stds[stds == 0] = 1.0
    z_matrix = (expr_matrix - means) / stds

    n_genes = len(merged)
    track_colors = [C_A if r["track"] == "A" else C_B for _, r in merged.iterrows()]
    gene_names = [label(neural, r[gid_col]) for _, r in merged.iterrows()]

    # Track sidebar
    ax_track.imshow([[1] for _ in range(n_genes)], aspect="auto", cmap="binary", vmin=0, vmax=1)
    for i, color in enumerate(track_colors):
        ax_track.add_patch(plt.Rectangle((-0.5, i - 0.5), 1, 1, color=color, ec="none"))
    ax_track.set_xticks([])
    ax_track.set_yticks([])
    ax_track.set_xlim(-0.5, 0.5)
    ax_track.set_ylim(n_genes - 0.5, -0.5)
    ax_track.set_ylabel("Track", fontsize=7.5, fontweight="bold")
    ax_track.spines[:].set_visible(False)

    # Main heatmap
    im = ax_heat.imshow(z_matrix, aspect="auto", cmap=plt.cm.coolwarm, vmin=-2.0, vmax=2.0, interpolation="nearest")
    ax_heat.set_xticks(range(len(time_labels)))
    ax_heat.set_xticklabels(time_labels, rotation=35, ha="right", fontsize=6.5)
    ax_heat.set_yticks(range(n_genes))
    ax_heat.set_yticklabels(gene_names, fontsize=7)
    ax_heat.axhline(4.5, color="#333333", lw=0.8, ls="--")
    ax_heat.set_xlabel("Regeneration timepoint", fontsize=7.5)
    ax_heat.set_ylabel("Candidate", fontsize=7.5)
    ax_heat.set_title("Temporal induction dynamics", fontsize=8, pad=4)
    ax_heat.spines[:].set_visible(False)
    panel_tag(ax_heat, "a")

    # Heatmap colorbar
    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("Relative expression (z-score)", fontsize=6.5)
    cbar.ax.tick_params(labelsize=6)

    # --- Panel b: Relative Fold-Induction Profiles ---
    time_x = np.arange(len(time_cols))
    # Compute fold change relative to 0d intact: Et / E0d
    fc_matrix = expr_matrix / (expr_matrix[:, [0]] + 1e-6)

    # Highlight exemplary early vs late factors
    early_genes = ["ascl-2", "ptf-4", "dd18719"]
    late_genes = ["dd22163", "Zeb-1", "dd18972"]

    for i, (_, r) in enumerate(merged.iterrows()):
        nm = gene_names[i]
        fc = fc_matrix[i]
        if nm in early_genes:
            ax_traj.plot(time_x, fc, "o-", color="#B04A3E", lw=1.5, markersize=3.5,
                        label=f"{nm} (early)", alpha=0.9, zorder=4)
        elif nm in late_genes:
            ax_traj.plot(time_x, fc, "s-", color="#2B4C6F", lw=1.5, markersize=3.5,
                        label=f"{nm} (late)", alpha=0.9, zorder=4)
        else:
            ax_traj.plot(time_x, fc, "-", color="#CCCCCC", lw=0.8, alpha=0.5, zorder=2)

    ax_traj.axhline(1.0, color="#888888", lw=0.8, ls=":", zorder=1)
    ax_traj.set_xticks(time_x)
    ax_traj.set_xticklabels(time_labels, rotation=35, ha="right", fontsize=6.5)
    ax_traj.set_xlabel("Regeneration timepoint", fontsize=7.5)
    ax_traj.set_ylabel("Fold induction ($E_t / E_{0d}$)", fontsize=7.5)
    ax_traj.set_title("Regeneration trajectory kinetics", fontsize=8, pad=4)
    ax_traj.set_ylim(0.2, 2.6)
    ax_traj.legend(loc="upper left", frameon=False, fontsize=5.8, ncol=2)
    ax_traj.spines["top"].set_visible(False)
    ax_traj.spines["right"].set_visible(False)
    panel_tag(ax_traj, "b")

    fig.tight_layout()
    save(fig, "36_regeneration_temporal_dynamics")
    print("  wrote 36_regeneration_temporal_dynamics")

if __name__ == "__main__":
    build()
