#!/usr/bin/env python
"""Visualize Dirichlet-centered (k=40) prioritization results — enhanced.

Generates 5 publication-quality figures from Dirichlet-centered outputs:

  1. trackA_top5          - Track A top-5 with domain, RNAi phenotype, score details
  2. trackB_top5          - Track B top-5 with domain, integrated score, score details
  3. scatter_fixed_vs_dirichlet - fixed vs Dirichlet composite with shift annotations
  4. combined_dual_track  - both tracks with composite, base score, domain annotations
  3. score_shift          - Dirichlet composite - fixed composite per candidate

Outputs: projects/NeuralTF/figures/fig_centered_*.png

Usage:
    python projects/NeuralTF/scripts/visualize_centered.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import spearmanr, pearsonr

REPO = Path(__file__).resolve().parents[3]
IN_DIR = REPO / "projects" / "NeuralTF" / "results"
OUT_DIR = REPO / "projects" / "NeuralTF" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Style — Nature journal standards
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 10.5,
    "axes.labelsize": 9.5,
    "axes.titleweight": "bold",
    "axes.labelweight": "bold",
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.12,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.spines.left": True,
    "axes.spines.bottom": True,
    "axes.edgecolor": "#333333",
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "text.color": "#333333",
    "legend.frameon": True,
    "legend.framealpha": 0.92,
    "legend.edgecolor": "0.8",
})

# Okabe-Ito colorblind-safe palette
C_A = "#D55E00"   # Track A (RNAi-validated)
C_B = "#0072B2"   # Track B (novel)
C_GRAY = "#999999"
C_GRID = "#E0E0E0"
C_BLACK = "#111111"
C_GREEN = "#009E73"
C_SKY = "#56B4E9"
C_ORANGE = "#E69F00"

DOMAIN_COLORS = {
    "bHLH": "#0072B2",
    "Homeobox": "#009E73",
    "Znf": "#CC79A7",
    "fork_head": "#E69F00",
    "T-box": "#D55E00",
    "Ets": "#56B4E9",
    "none": "#BBBBBB",
}


def _style_ax(ax):
    ax.xaxis.grid(True, color=C_GRID, linewidth=0.5, zorder=0)
    ax.yaxis.grid(True, color=C_GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(colors=C_BLACK)
    for spine in ax.spines.values():
        spine.set_color(C_BLACK)
    return ax


def _add_panel_letter(ax, letter, x=-0.08, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top", ha="left", color=C_BLACK)


def _domain_group(domains: str) -> str:
    if not isinstance(domains, str) or not domains.strip():
        return "none"
    d = domains.lower()
    if "bhlh" in d: return "bHLH"
    if "homeobox" in d: return "Homeobox"
    if "znf" in d: return "Znf"
    if "fork_head" in d: return "fork_head"
    if "t-box" in d: return "T-box"
    if "ets" in d: return "Ets"
    return "none"


def _add_panel_letter(ax, letter, x=-0.08, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top", ha="left", color=C_BLACK)


# ---------------------------------------------------------------------------
# Figure builders
# ---------------------------------------------------------------------------
def fig_trackA_top5() -> None:
    csv_path = IN_DIR / "dirichlet_top10_prioritized.csv"
    if not csv_path.exists():
        return
    top = pd.read_csv(csv_path)
    a = top[top["track"] == "A"].copy()
    if a.empty:
        return
    a = a.sort_values("dirichlet_median_score", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8.5, 4))
    _add_panel_letter(ax, "A")
    y = np.arange(len(a))
    dom_groups = [_domain_group(d) for d in a["interpro_domains"]]
    colors = [DOMAIN_COLORS.get(g, C_GRAY) for g in dom_groups]

    bars = ax.barh(y, a["composite_score"], color=colors, edgecolor="white",
                   height=0.6, linewidth=0.5, zorder=3)
    ax.set_yticks(y)
    # Enhanced labels: gene name + dirichlet median + domain
    labels = [f"{r['gene_name']}  (median={r['dirichlet_median_score']:.3f}, {_domain_group(r['interpro_domains'])})" 
              for _, r in a.iterrows()]
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlim(0, max(1.05, a["composite_score"].max() * 1.15))
    ax.set_xlabel("Composite score (Dirichlet median + bonuses)", fontsize=10)
    ax.set_title("A · Track A — RNAi-validated (Dirichlet-centered, k=40)",
                 fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    _style_ax(ax)

    # Add score values on bars
    for bar, score in zip(bars, a["composite_score"]):
        ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2, 
                f"{score:.3f}", va='center', fontsize=7, fontweight='bold')

    # Legend
    handles = [Patch(facecolor=DOMAIN_COLORS[d], edgecolor="white",
                     label=d) for d in sorted(set(dom_groups)) if d != "none"]
    ax.legend(handles=handles, loc="lower right", fontsize=8, title="Domain group", framealpha=0.9)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_centered_trackA_top5.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_centered_trackA_top5.png")


def fig_trackB_top5() -> None:
    csv_path = IN_DIR / "dirichlet_top10_prioritized.csv"
    if not csv_path.exists():
        return
    top = pd.read_csv(csv_path)
    b = top[top["track"] == "B"].copy()
    if b.empty:
        return
    b = b.sort_values("dirichlet_median_score", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8.5, 4))
    _add_panel_letter(ax, "B")
    y = np.arange(len(b))
    dom_groups = [_domain_group(d) for d in b["interpro_domains"]]
    colors = [DOMAIN_COLORS.get(g, C_GRAY) for g in dom_groups]

    bars = ax.barh(y, b["composite_score"], color=colors, edgecolor="white",
                   height=0.6, linewidth=0.5, zorder=3)
    ax.set_yticks(y)
    labels = [f"{r['gene_name']}  (median={r['dirichlet_median_score']:.3f}, {_domain_group(r['interpro_domains'])})" 
              for _, r in b.iterrows()]
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlim(0, max(1.05, b["composite_score"].max() * 1.15))
    ax.set_xlabel("Composite score (Dirichlet median + bonuses)", fontsize=10)
    ax.set_title("B · Track B — Novel candidates (Dirichlet-centered, k=40)",
                 fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    _style_ax(ax)

    for bar, score in zip(bars, b["composite_score"]):
        ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2, 
                f"{score:.3f}", va='center', fontsize=7, fontweight='bold')

    handles = [Patch(facecolor=DOMAIN_COLORS[d], edgecolor="white",
                     label=d) for d in sorted(set(dom_groups)) if d != "none"]
    ax.legend(handles=handles, loc="lower right", fontsize=8, title="Domain group", framealpha=0.9)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_centered_trackB_top5.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_centered_trackB_top5.png")


def fig_scatter_fixed_vs_dirichlet() -> None:
    centered_csv = IN_DIR / "dirichlet_top10_prioritized.csv"
    fixed_csv = IN_DIR / "top10_neural_tfs_prioritized.csv"

    if not centered_csv.exists() or not fixed_csv.exists():
        return

    centered = pd.read_csv(centered_csv)
    fixed = pd.read_csv(fixed_csv)

    merged = pd.merge(fixed[["gene_id_v6", "composite_score"]],
                      centered[["gene_id_v6", "composite_score", "dirichlet_median_score", "track"]],
                      on="gene_id_v6", suffixes=("_fixed", "_dirichlet"))

    merged["shift"] = merged["composite_score_dirichlet"] - merged["composite_score_fixed"]
    rho, p_val = spearmanr(merged["composite_score_fixed"], merged["composite_score_dirichlet"])
    r_pearson, p_pearson = pearsonr(merged["composite_score_fixed"], merged["composite_score_dirichlet"])

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    _add_panel_letter(ax, "C")

    for _, row in merged.iterrows():
        color = C_A if row["track"] == "A" else C_B
        marker = "o" if row["track"] == "A" else "^"
        ax.scatter(row["composite_score_fixed"], row["composite_score_dirichlet"],
                   color=color, s=100, marker=marker, edgecolors="white",
                   linewidth=1, zorder=3)
        ax.annotate(row["gene_id_v6"].replace("dd_Smed_v6_", "dd"),
                    (row["composite_score_fixed"], row["composite_score_dirichlet"]),
                    textcoords="offset points", xytext=(6, 3),
                    fontsize=7, color=C_BLACK, fontweight='bold')

    lims = [0, 1.0]
    ax.plot(lims, lims, "--", color=C_GRAY, linewidth=1.5, zorder=1,
            label="y = x (no shift)")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Fixed-weight composite score", fontsize=10)
    ax.set_ylabel("Dirichlet-centered composite score", fontsize=10)
    ax.set_title(f"C · Fixed vs Dirichlet-centered (ρ={rho:.3f}, r={r_pearson:.3f}, p={p_val:.2e})",
                 fontsize=11, fontweight="bold")
    
    # Legend with correlation info
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_A, markersize=10, label="Track A (RNAi)"),
        Line2D([0], [0], marker="^", color="w", markerfacecolor=C_B, markersize=10, label="Track B (novel)"),
        Line2D([0], [0], color=C_GRAY, linewidth=1.5, linestyle="--", label="y = x (no shift)"),
        Line2D([0], [0], color="white", label=f"Spearman ρ = {rho:.3f}"),
        Line2D([0], [0], color="white", label=f"Pearson r = {r_pearson:.3f}"),
    ]
    ax.legend(handles=handles, fontsize=7, loc="upper left", framealpha=0.9)
    _style_ax(ax)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_centered_scatter_fixed_vs_dirichlet.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_centered_scatter_fixed_vs_dirichlet.png")


def fig_combined_dual_track() -> None:
    csv_path = IN_DIR / "dirichlet_top10_prioritized.csv"
    if not csv_path.exists():
        return
    top = pd.read_csv(csv_path)
    a = top[top["track"] == "A"].copy()
    b = top[top["track"] == "B"].copy()

    fig, ax = plt.subplots(figsize=(9, 5.5))
    _add_panel_letter(ax, "D")

    for track_df, color, label, marker in [(a, C_A, "Track A (RNAi)", "o"),
                                             (b, C_B, "Track B (novel)", "^")]:
        if track_df.empty:
            continue
        track_df = track_df.sort_values("composite_score", ascending=False).reset_index(drop=True)
        x = np.arange(len(track_df))
        dom_groups = [_domain_group(d) for d in track_df["interpro_domains"]]
        colors = [DOMAIN_COLORS.get(g, C_GRAY) for g in dom_groups]
        bars = ax.bar(x, track_df["composite_score"], color=colors, edgecolor="white",
                      width=0.6, linewidth=0.5, zorder=3, label=label)
        for i, (_, row) in enumerate(track_df.iterrows()):
            # Show gene name + composite score
            ax.text(i, row["composite_score"] + 0.012, row["gene_name"],
                    ha="center", va="bottom", fontsize=7, rotation=0)
            # Add composite score value on bar
            ax.text(i, row["composite_score"]/2, f"{row['composite_score']:.3f}",
                    ha="center", va="center", fontsize=7, color="white", fontweight='bold')

    ax.set_xticks([])
    ax.set_ylabel("Composite score (Dirichlet median + bonuses)", fontsize=10)
    ax.set_ylim(0, max(top["composite_score"].max() * 1.25, 1.0))
    ax.set_title("D · Dirichlet-centered dual track — Top-10", 
                 fontsize=11, fontweight="bold")

    handles1 = [Line2D([0], [0], marker="o", color="w", markerfacecolor=C_A,
                       markersize=10, label="Track A (RNAi)"),
                Line2D([0], [0], marker="^", color="w", markerfacecolor=C_B,
                       markersize=10, label="Track B (novel)")]
    handles2 = [Patch(facecolor=DOMAIN_COLORS[d], edgecolor="white", label=d)
                for d in sorted(set(_domain_group(d) for d in top["interpro_domains"])) if d != "none"]
    leg1 = ax.legend(handles=handles1, loc="upper left", fontsize=8, framealpha=0.9)
    leg2 = ax.legend(handles=handles2, loc="upper right", fontsize=8, title="Domain group", framealpha=0.9)
    ax.add_artist(leg1)
    _style_ax(ax)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_centered_combined_dual_track.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_centered_combined_dual_track.png")


def fig_score_shift() -> None:
    centered_csv = IN_DIR / "dirichlet_top10_prioritized.csv"
    fixed_csv = IN_DIR / "top10_neural_tfs_prioritized.csv"

    if not centered_csv.exists() or not fixed_csv.exists():
        return

    centered = pd.read_csv(centered_csv)
    fixed = pd.read_csv(fixed_csv)

    merged = pd.merge(fixed[["gene_id_v6", "composite_score"]],
                      centered[["gene_id_v6", "composite_score", "track"]],
                      on="gene_id_v6", suffixes=("_fixed", "_dirichlet"))

    merged["shift"] = merged["composite_score_dirichlet"] - merged["composite_score_fixed"]
    merged = merged.sort_values("shift", ascending=False).reset_index(drop=True)
    
    mean_shift = merged["shift"].mean()
    pos_shifts = (merged["shift"] > 0).sum()
    neg_shifts = (merged["shift"] < 0).sum()

    fig, ax = plt.subplots(figsize=(9, 5.5))
    _add_panel_letter(ax, "E")

    colors = [C_A if t == "A" else C_B for t in merged["track"]]
    x = np.arange(len(merged))
    bars = ax.bar(x, merged["shift"], color=colors, edgecolor="white",
                  width=0.6, linewidth=0.5, zorder=3)
    ax.axhline(0, color=C_BLACK, linewidth=1.0, linestyle="--")
    ax.axhline(mean_shift, color=C_ORANGE, linewidth=1.5, linestyle=":", 
               label=f"Mean shift = {mean_shift:+.4f}")
    ax.set_xticks(x)
    ax.set_xticklabels(merged["gene_id_v6"].str.replace("dd_Smed_v6_", "dd"),
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Dirichlet composite − Fixed composite", fontsize=10)
    ax.set_title(f"E · Score shift (Dirichlet-centered − fixed) | Mean: {mean_shift:+.4f} | ↑{pos_shifts} ↓{neg_shifts}",
                 fontsize=11, fontweight="bold")
    ax.legend(handles=[
        Patch(facecolor=C_A, label="Track A (RNAi)"),
        Patch(facecolor=C_B, label="Track B (novel)"),
        Line2D([0], [0], color=C_ORANGE, linewidth=1.5, linestyle=":", label=f"Mean = {mean_shift:+.4f}"),
    ], loc="upper right", fontsize=8, framealpha=0.9)
    _style_ax(ax)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_centered_score_shift.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_centered_score_shift.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
IN_DIR = Path(__file__).resolve().parents[3] / "projects" / "NeuralTF" / "results"
OUT_DIR = Path(__file__).resolve().parents[3] / "projects" / "NeuralTF" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Style — Nature journal standards
# ---------------------------------------------------------------------------
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 10.5,
    "axes.labelsize": 9.5,
    "axes.titleweight": "bold",
    "axes.labelweight": "bold",
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.12,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.spines.left": True,
    "axes.spines.bottom": True,
    "axes.edgecolor": "#333333",
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "text.color": "#333333",
    "legend.frameon": True,
    "legend.framealpha": 0.92,
    "legend.edgecolor": "0.8",
})

# Okabe-Ito colorblind-safe palette
C_A = "#D55E00"
C_B = "#0072B2"
C_GRAY = "#999999"
C_GRID = "#E0E0E0"
C_BLACK = "#111111"
C_GREEN = "#009E73"
C_SKY = "#56B4E9"
C_ORANGE = "#E69F00"

DOMAIN_COLORS = {
    "bHLH": "#0072B2",
    "Homeobox": "#009E73",
    "Znf": "#CC79A7",
    "fork_head": "#E69F00",
    "T-box": "#D55E00",
    "Ets": "#56B4E9",
    "none": "#BBBBBB",
}


def _style_ax(ax):
    ax.xaxis.grid(True, color=C_GRID, linewidth=0.5, zorder=0)
    ax.yaxis.grid(True, color=C_GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(colors=C_BLACK)
    for spine in ax.spines.values():
        spine.set_color(C_BLACK)
    return ax


def _add_panel_letter(ax, letter, x=-0.08, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top", ha="left", color=C_BLACK)


def _domain_group(domains: str) -> str:
    if not isinstance(domains, str) or not domains.strip():
        return "none"
    d = domains.lower()
    if "bhlh" in d: return "bHLH"
    if "homeobox" in d: return "Homeobox"
    if "znf" in d: return "Znf"
    if "fork_head" in d: return "fork_head"
    if "t-box" in d: return "T-box"
    if "ets" in d: return "Ets"
    return "none"


def _add_panel_letter(ax, letter, x=-0.08, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top", ha="left", color=C_BLACK)


def main() -> int:
    print("== Visualize Dirichlet-centered (k=40) prioritization ==")
    fig_trackA_top5()
    fig_trackB_top5()
    fig_scatter_fixed_vs_dirichlet()
    fig_combined_dual_track()
    fig_score_shift()
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())