"""Figure 6 — Prioritized candidate atlas.

The biological synthesis figure: compact visual encoding of all evidence
for the final Top 10 candidates.

Panels:
  A  Candidate information card (integrated visual encoding)
  B  Multi-method score comparison (lollipop plot)
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D


def _load_data():
    neural_df = load_rank_neural()
    fixed = load_top10_fixed()
    centered = load_top10_centered()
    uniform = load_top10_uniform()
    return neural_df, fixed, centered, uniform


def _get_candidate_info(gene_id, neural_df, fixed, centered, uniform):
    """Gather all information for a single candidate."""
    info = {"gene_id": gene_id}

    # From fixed top10
    f_row = fixed[fixed["gene_id"] == gene_id]
    if len(f_row) > 0:
        info["track"] = f_row.iloc[0].get("track", "")
        info["gene_name"] = f_row.iloc[0].get("gene_name", "")
        info["composite_fixed"] = f_row.iloc[0].get("composite_score", np.nan)
        info["proof_status"] = f_row.iloc[0].get("proof_status", "")
        info["domains"] = f_row.iloc[0].get("interpro_domains",
                           f_row.iloc[0].get("domains_all", ""))
        info["human_ortholog"] = f_row.iloc[0].get("human_ortholog", "")
        info["rnai_phenotype"] = f_row.iloc[0].get("rnai_phenotype_notes", "")

    # From centered
    c_row = centered[centered["gene_id"] == gene_id]
    if len(c_row) > 0:
        info["composite_centered"] = c_row.iloc[0].get("composite_score", np.nan)
        info["dirichlet_median"] = c_row.iloc[0].get("dirichlet_median_score", np.nan)

    # From uniform
    u_row = uniform[uniform["gene_id"] == gene_id]
    if len(u_row) > 0:
        info["composite_uniform"] = u_row.iloc[0].get("composite_score", np.nan)
        info["uniform_median"] = u_row.iloc[0].get("uniform_median_score", np.nan)

    # From neural_df (evidence streams)
    n_row = neural_df[neural_df["gene_id"] == gene_id]
    if len(n_row) > 0:
        for col in STREAM_COLS:
            info[col] = n_row.iloc[0].get(col, np.nan)
        info["integrated_score"] = n_row.iloc[0].get("integrated_score", np.nan)

    return info


STREAM_COLS = ["expression", "specificity", "reproducibility",
               "rnai", "correlation", "neural_enriched", "neural_specificity"]


def fig6a_candidate_cards(neural_df, fixed, centered, uniform, ax):
    """Compact candidate information cards for Top 10."""
    ax.axis("off")

    # Order candidates
    ordered = fixed.sort_values("composite_score", ascending=False)
    candidates = []
    for _, row in ordered.iterrows():
        info = _get_candidate_info(row["gene_id"], neural_df, fixed, centered, uniform)
        candidates.append(info)

    # Layout: 2 columns × 5 rows
    n_cols = 2
    n_rows = 5

    for idx, info in enumerate(candidates):
        col = idx % n_cols
        row = idx // n_rows

        x_start = col * 0.5
        y_start = 1.0 - (idx // n_cols) * 0.18 - 0.02

        track = info.get("track", "")
        track_color = C_TRACK_A if track == "A" else C_TRACK_B
        track_symbol = "●" if track == "A" else "■"

        # Gene ID and name
        name = info.get("gene_name", "")
        gid = info.get("gene_id", "")
        label = f"{gid}" if not name else f"{name} ({gid})"

        ax.text(x_start, y_start, f"{track_symbol} {label}",
                fontsize=6, fontweight="bold", color=track_color,
                transform=ax.transAxes, va="top")

        # Score bar
        composite = info.get("composite_fixed", 0)
        bar_x = x_start + 0.02
        bar_w = composite * 0.35  # Scale to fit
        ax.add_patch(FancyBboxPatch(
            (bar_x, y_start - 0.035), bar_w, 0.02,
            boxstyle="round,pad=0.002", facecolor=track_color, alpha=0.6,
            transform=ax.transAxes))

        ax.text(bar_x + bar_w + 0.01, y_start - 0.025, f"{composite:.3f}",
                fontsize=5, va="center", transform=ax.transAxes, color="#333333")

        # Proof status
        proof = info.get("proof_status", "")
        if "validated" in str(proof).lower() or "fstf" in str(proof).lower():
            proof_label = "[validated]"
            proof_color = "#0072B2"
        else:
            proof_label = "[novel]"
            proof_color = "#E69F00"
        ax.text(x_start + 0.42, y_start, proof_label,
                fontsize=5, color=proof_color, transform=ax.transAxes, va="top")

        # Evidence streams as tiny colored dots
        dot_x = x_start + 0.02
        dot_y = y_start - 0.06
        for j, stream in enumerate(STREAM_COLS):
            val = info.get(stream, np.nan)
            if pd.notna(val) and val > 0:
                alpha = 0.3 + 0.7 * val
                ax.plot(dot_x + j * 0.025, dot_y, "o",
                        color=STREAM_COLORS[stream], markersize=2.5,
                        alpha=alpha, transform=ax.transAxes)
            else:
                ax.plot(dot_x + j * 0.025, dot_y, "o",
                        color="#DDDDDD", markersize=2.5,
                        alpha=0.5, transform=ax.transAxes)

        # Domain annotation
        domains = str(info.get("domains", ""))
        if domains and domains != "nan":
            # Show first domain only
            first_domain = domains.split(";")[0].strip()[:25]
            ax.text(x_start + 0.02, y_start - 0.09, f"Domain: {first_domain}",
                    fontsize=4.5, color="#666666", transform=ax.transAxes, va="top",
                    fontstyle="italic")

        # Ortholog
        ortholog = str(info.get("human_ortholog", ""))
        if ortholog and ortholog != "nan":
            ax.text(x_start + 0.25, y_start - 0.09, f"Hsa: {ortholog}",
                    fontsize=4.5, color="#666666", transform=ax.transAxes, va="top")

    # Stream legend at bottom
    for j, stream in enumerate(STREAM_COLS):
        ax.plot(0.02 + j * 0.07, 0.02, "o", color=STREAM_COLORS[stream],
                markersize=4, transform=ax.transAxes)
        ax.text(0.04 + j * 0.07, 0.02, STREAM_LABELS[stream][:4],
                fontsize=4.5, va="center", transform=ax.transAxes, color="#666666")

    ax.set_title("Prioritized candidate atlas", fontsize=8, fontweight="bold",
                 pad=6, loc="left")


def fig6b_multi_method(fixed, centered, uniform, ax):
    """Lollipop plot comparing scores across three methods for Top 10."""
    # Build comparison table
    candidates = []
    for _, row in fixed.iterrows():
        gid = row["gene_id"]
        track = row.get("track", "")
        f_score = row.get("composite_score", np.nan)

        c_row = centered[centered["gene_id"] == gid]
        c_score = c_row.iloc[0].get("composite_score", np.nan) if len(c_row) > 0 else np.nan

        u_row = uniform[uniform["gene_id"] == gid]
        u_score = u_row.iloc[0].get("composite_score", np.nan) if len(u_row) > 0 else np.nan

        candidates.append({
            "gene_id": gid, "track": track,
            "fixed": f_score, "centered": c_score, "uniform": u_score,
            "label": gene_label(fixed, gid),
        })

    df = pd.DataFrame(candidates)
    df = df.sort_values("fixed", ascending=False)

    y = np.arange(len(df))[::-1]
    methods = ["fixed", "centered", "uniform"]
    method_colors = {"fixed": C_FIXED, "centered": C_CENTERED, "uniform": C_UNIFORM}
    offsets = [-0.2, 0, 0.2]

    for method, offset in zip(methods, offsets):
        vals = df[method].values
        valid = ~np.isnan(vals)
        ax.hlines(y[valid] + offset, 0, vals[valid],
                  color=method_colors[method], linewidth=0.8, alpha=0.7)
        ax.scatter(vals[valid], y[valid] + offset, color=method_colors[method],
                   s=8, zorder=5, edgecolors="white", linewidth=0.3)

    # Y-axis labels
    y_labels = []
    for _, row in df.iterrows():
        track = row["track"]
        prefix = "●" if track == "A" else "■"
        y_labels.append(f"{prefix} {row['label']}")
    ax.set_yticks(y)
    ax.set_yticklabels(y_labels, fontsize=5.5)

    # Color y-tick labels
    for i, (_, row) in enumerate(df.iterrows()):
        color = C_TRACK_A if row["track"] == "A" else C_TRACK_B
        ax.get_yticklabels()[i].set_color(color)

    ax.set_xlabel("Composite score")
    style_ax(ax, title="")
    ax.set_title("Cross-method comparison", fontsize=8, fontweight="bold",
                 pad=6, loc="left")

    # Legend
    legend_elements = [
        Line2D([0], [0], color=C_FIXED, marker="o", label="Fixed", linewidth=0.8, markersize=4),
        Line2D([0], [0], color=C_CENTERED, marker="o", label="Centered", linewidth=0.8, markersize=4),
        Line2D([0], [0], color=C_UNIFORM, marker="o", label="Uniform", linewidth=0.8, markersize=4),
    ]
    ax.legend(handles=legend_elements, frameon=False, fontsize=5, loc="lower right")


def build():
    """Generate Figure 6 panels A–B."""
    neural_df, fixed, centered, uniform = _load_data()

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL, inch(80)),
                             gridspec_kw={"width_ratios": [1.5, 1.0]})

    fig6a_candidate_cards(neural_df, fixed, centered, uniform, axes[0])
    fig6b_multi_method(fixed, centered, uniform, axes[1])

    panel_letter(axes[0], "A")
    panel_letter(axes[1], "B")

    fig.tight_layout(w_pad=1.5)
    savefig(fig, "Fig6_prioritized_candidate_atlas")
    return fig


if __name__ == "__main__":
    build()
