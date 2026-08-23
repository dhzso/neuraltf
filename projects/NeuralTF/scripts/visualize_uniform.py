#!/usr/bin/env python
"""Visualize Dirichlet-uniform (alpha=1) prioritization results.

Generates 7 publication-quality figures from Dirichlet-uniform outputs:

  1. trackA_top5          - Track A top-5 (RNAi-validated) with domain colors
  2. trackB_top5          - Track B top-5 (novel) with domain colors
  3. scatter_fixed_vs_uniform - fixed-weight vs uniform median
  4. scatter_centered_vs_uniform - centered vs uniform median
  5. combined_dual_track  - both tracks with composite overlay
  6. score_shift          - Uniform composite - fixed composite per candidate
  7. three_way_comparison - grouped bars: fixed / centered / uniform

Outputs: projects/NeuralTF/figures/fig_uniform_*.png

Usage:
    python projects/NeuralTF/scripts/visualize_uniform.py
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
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.12,
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

C_FIXED = "#666666"
C_CENTERED = "#E69F00"
C_UNIFORM = "#009E73"


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
    if "bhlh" in d:
        return "bHLH"
    if "homeobox" in d:
        return "Homeobox"
    if "znf" in d:
        return "Znf"
    if "fork_head" in d:
        return "fork_head"
    if "t-box" in d:
        return "T-box"
    if "ets" in d:
        return "Ets"
    return "none"


# ---------------------------------------------------------------------------
# Figure builders
# ---------------------------------------------------------------------------
def fig_trackA_top5() -> None:
    csv_path = IN_DIR / "dirichlet_uniform_top10.csv"
    if not csv_path.exists():
        return
    top = pd.read_csv(csv_path)
    a = top[top["track"] == "A"].copy()
    if a.empty:
        return
    a = a.sort_values("uniform_median_score", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(7.8, 3.6))
    _add_panel_letter(ax, "A")
    y = np.arange(len(a))
    dom_groups = [_domain_group(d) for d in a["interpro_domains"]]
    colors = [DOMAIN_COLORS.get(g, C_GRAY) for g in dom_groups]

    ax.barh(y, a["composite_score"], color=colors, edgecolor="white",
            height=0.6, linewidth=0.5, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(a["gene_name"], fontsize=9.5)
    ax.set_xlim(0, max(1.05, a["composite_score"].max() * 1.15))
    ax.set_xlabel("Composite score (Uniform median + bonuses)", fontsize=10)
    ax.set_title("A · Track A — RNAi-validated (Dirichlet-uniform)",
                 fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    _style_ax(ax)

    handles = [Patch(facecolor=DOMAIN_COLORS[d], edgecolor="white",
                     label=d) for d in sorted(set(dom_groups)) if d != "none"]
    ax.legend(handles=handles, loc="lower right", fontsize=8, title="Domain group")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_uniform_trackA_top5.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_uniform_trackA_top5.png")


def fig_trackB_top5() -> None:
    csv_path = IN_DIR / "dirichlet_uniform_top10.csv"
    if not csv_path.exists():
        return
    top = pd.read_csv(csv_path)
    b = top[top["track"] == "B"].copy()
    if b.empty:
        return
    b = b.sort_values("uniform_median_score", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(7.8, 3.6))
    _add_panel_letter(ax, "B")
    y = np.arange(len(b))
    dom_groups = [_domain_group(d) for d in b["interpro_domains"]]
    colors = [DOMAIN_COLORS.get(g, C_GRAY) for g in dom_groups]

    ax.barh(y, b["composite_score"], color=colors, edgecolor="white",
            height=0.6, linewidth=0.5, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(b["gene_name"], fontsize=9.5)
    ax.set_xlim(0, max(1.05, b["composite_score"].max() * 1.15))
    ax.set_xlabel("Composite score (Uniform median + bonuses)", fontsize=10)
    ax.set_title("B · Track B — Novel candidates (Dirichlet-uniform)",
                 fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    _style_ax(ax)

    handles = [Patch(facecolor=DOMAIN_COLORS[d], edgecolor="white",
                     label=d) for d in sorted(set(dom_groups)) if d != "none"]
    ax.legend(handles=handles, loc="lower right", fontsize=8, title="Domain group")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_uniform_trackB_top5.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_uniform_trackB_top5.png")


def fig_scatter_fixed_vs_uniform() -> None:
    uniform_csv = IN_DIR / "dirichlet_uniform_top10.csv"
    fixed_csv = IN_DIR / "top10_neural_tfs_prioritized.csv"

    if not uniform_csv.exists() or not fixed_csv.exists():
        return

    uniform = pd.read_csv(uniform_csv)
    fixed = pd.read_csv(fixed_csv)

    merged = pd.merge(fixed[["gene_id_v6", "composite_score"]],
                      uniform[["gene_id_v6", "composite_score", "uniform_median_score", "track"]],
                      on="gene_id_v6", suffixes=("_fixed", "_uniform"))

    fig, ax = plt.subplots(figsize=(6, 5))
    _add_panel_letter(ax, "C")

    for _, row in merged.iterrows():
        color = C_A if row["track"] == "A" else C_B
        marker = "o" if row["track"] == "A" else "^"
        ax.scatter(row["composite_score_fixed"], row["composite_score_uniform"],
                   color=color, s=80, marker=marker, edgecolors="white",
                   linewidth=0.8, zorder=3)
        ax.annotate(row["gene_id_v6"].replace("dd_Smed_v6_", "dd"),
                    (row["composite_score_fixed"], row["composite_score_uniform"]),
                    textcoords="offset points", xytext=(5, 3),
                    fontsize=7, color=C_BLACK)

    lims = [0, 1.0]
    ax.plot(lims, lims, "--", color=C_GRAY, linewidth=1.0, zorder=1,
            label="y = x (no shift)")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Fixed-weight composite score", fontsize=10)
    ax.set_ylabel("Dirichlet-uniform composite score", fontsize=10)
    ax.set_title("C · Fixed-weight vs Dirichlet-uniform composite",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    _style_ax(ax)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_uniform_scatter_fixed_vs_uniform.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_uniform_scatter_fixed_vs_uniform.png")


def fig_scatter_centered_vs_uniform() -> None:
    uniform_csv = IN_DIR / "dirichlet_uniform_top10.csv"
    centered_csv = IN_DIR / "dirichlet_top10_prioritized.csv"

    if not uniform_csv.exists() or not centered_csv.exists():
        return

    uniform = pd.read_csv(uniform_csv)
    centered = pd.read_csv(centered_csv)

    merged = pd.merge(centered[["gene_id_v6", "composite_score", "track"]],
                      uniform[["gene_id_v6", "composite_score", "uniform_median_score"]],
                      on="gene_id_v6", suffixes=("_centered", "_uniform"))

    fig, ax = plt.subplots(figsize=(6, 5))
    _add_panel_letter(ax, "D")

    for _, row in merged.iterrows():
        color = C_A if row["track"] == "A" else C_B
        marker = "o" if row["track"] == "A" else "^"
        ax.scatter(row["composite_score_centered"], row["composite_score_uniform"],
                   color=color, s=80, marker=marker, edgecolors="white",
                   linewidth=0.8, zorder=3)
        ax.annotate(row["gene_id_v6"].replace("dd_Smed_v6_", "dd"),
                    (row["composite_score_centered"], row["composite_score_uniform"]),
                    textcoords="offset points", xytext=(5, 3),
                    fontsize=7, color=C_BLACK)

    lims = [0, 1.0]
    ax.plot(lims, lims, "--", color=C_GRAY, linewidth=1.0, zorder=1,
            label="y = x (no shift)")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Dirichlet-centered composite score", fontsize=10)
    ax.set_ylabel("Dirichlet-uniform composite score", fontsize=10)
    ax.set_title("D · Dirichlet-centered vs Dirichlet-uniform composite",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    _style_ax(ax)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_uniform_scatter_centered_vs_uniform.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_uniform_scatter_centered_vs_uniform.png")


def fig_combined_dual_track() -> None:
    csv_path = IN_DIR / "dirichlet_uniform_top10.csv"
    if not csv_path.exists():
        return
    top = pd.read_csv(csv_path)
    a = top[top["track"] == "A"].copy()
    b = top[top["track"] == "B"].copy()

    fig, ax = plt.subplots(figsize=(8, 5))
    _add_panel_letter(ax, "E")

    for track_df, color, label, marker in [(a, C_A, "Track A (RNAi)", "o"),
                                             (b, C_B, "Track B (novel)", "^")]:
        if track_df.empty:
            continue
        track_df = track_df.sort_values("composite_score", ascending=False).reset_index(drop=True)
        x = np.arange(len(track_df))
        dom_groups = [_domain_group(d) for d in track_df["interpro_domains"]]
        colors = [DOMAIN_COLORS.get(g, C_GRAY) for g in dom_groups]
        ax.bar(x, track_df["composite_score"], color=colors, edgecolor="white",
               width=0.6, linewidth=0.5, zorder=3, label=label)
        for i, (_, row) in enumerate(track_df.iterrows()):
            ax.text(i, row["composite_score"] + 0.015, row["gene_name"],
                    ha="center", va="bottom", fontsize=7, rotation=0)

    ax.set_xticks([])
    ax.set_ylabel("Composite score (Uniform median + bonuses)", fontsize=10)
    ax.set_ylim(0, max(top["composite_score"].max() * 1.2, 1.0))
    ax.set_title("E · Dirichlet-uniform dual track — Top-10",
                 fontsize=11, fontweight="bold")

    handles1 = [Line2D([0], [0], marker="o", color="w", markerfacecolor=C_A,
                       markersize=10, label="Track A (RNAi)"),
                Line2D([0], [0], marker="^", color="w", markerfacecolor=C_B,
                       markersize=10, label="Track B (novel)")]
    handles2 = [Patch(facecolor=DOMAIN_COLORS[d], edgecolor="white", label=d)
                for d in sorted(set(_domain_group(d) for d in top["interpro_domains"])) if d != "none"]
    leg1 = ax.legend(handles=handles1, loc="upper left", fontsize=8)
    leg2 = ax.legend(handles=handles2, loc="upper right", fontsize=8, title="Domain group")
    ax.add_artist(leg1)
    _style_ax(ax)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_uniform_combined_dual_track.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_uniform_combined_dual_track.png")


def fig_score_shift() -> None:
    uniform_csv = IN_DIR / "dirichlet_uniform_top10.csv"
    fixed_csv = IN_DIR / "top10_neural_tfs_prioritized.csv"

    if not uniform_csv.exists() or not fixed_csv.exists():
        return

    uniform = pd.read_csv(uniform_csv)
    fixed = pd.read_csv(fixed_csv)

    merged = pd.merge(fixed[["gene_id_v6", "composite_score"]],
                      uniform[["gene_id_v6", "composite_score", "track"]],
                      on="gene_id_v6", suffixes=("_fixed", "_uniform"))

    merged["shift"] = merged["composite_score_uniform"] - merged["composite_score_fixed"]
    merged = merged.sort_values("shift", ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    _add_panel_letter(ax, "F")

    colors = [C_A if t == "A" else C_B for t in merged["track"]]
    x = np.arange(len(merged))
    ax.bar(x, merged["shift"], color=colors, edgecolor="white",
           width=0.6, linewidth=0.5, zorder=3)
    ax.axhline(0, color=C_BLACK, linewidth=0.8, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(merged["gene_id_v6"].str.replace("dd_Smed_v6_", "dd"),
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Uniform composite − Fixed composite", fontsize=10)
    ax.set_title("F · Score shift (Dirichlet-uniform − fixed)",
                 fontsize=11, fontweight="bold")
    ax.legend(handles=[
        Patch(facecolor=C_A, label="Track A (RNAi)"),
        Patch(facecolor=C_B, label="Track B (novel)")
    ], loc="upper right", fontsize=8)
    _style_ax(ax)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_uniform_score_shift.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_uniform_score_shift.png")


def fig_three_way_comparison() -> None:
    """Grouped bar chart: fixed / centered / uniform for all top-10 candidates."""
    fixed_csv = IN_DIR / "top10_neural_tfs_prioritized.csv"
    centered_csv = IN_DIR / "dirichlet_top10_prioritized.csv"
    uniform_csv = IN_DIR / "dirichlet_uniform_top10.csv"

    if not all(p.exists() for p in [fixed_csv, centered_csv, uniform_csv]):
        return

    fixed = pd.read_csv(fixed_csv)
    centered = pd.read_csv(centered_csv)
    uniform = pd.read_csv(uniform_csv)

    # Merge all three on gene_id_v6
    merged = fixed[["gene_id_v6", "composite_score"]].rename(
        columns={"composite_score": "fixed"})
    merged = merged.merge(centered[["gene_id_v6", "composite_score"]].rename(
        columns={"composite_score": "centered"}), on="gene_id_v6")
    merged = merged.merge(uniform[["gene_id_v6", "composite_score", "track"]].rename(
        columns={"composite_score": "uniform"}), on="gene_id_v6")

    # Sort by fixed composite score
    merged = merged.sort_values("fixed", ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10, 5))
    _add_panel_letter(ax, "G")

    x = np.arange(len(merged))
    width = 0.25
    colors = [C_FIXED, C_CENTERED, C_UNIFORM]
    labels = ["Fixed-weight", "Dirichlet-centered (k=40)", "Dirichlet-uniform (α=1)"]

    for i, (col, color, label) in enumerate(zip(["fixed", "centered", "uniform"],
                                                   colors, labels)):
        ax.bar(x + (i - 1) * width, merged[col], width,
               color=color, edgecolor="white", linewidth=0.5,
               label=label, zorder=3)

    # Track color coding on x-axis
    for i, (_, row) in enumerate(merged.iterrows()):
        track_color = C_A if row["track"] == "A" else C_B
        ax.text(i, -0.03, row["gene_id_v6"].replace("dd_Smed_v6_", "dd"),
                ha="center", va="top", fontsize=7, color=track_color,
                transform=ax.get_xaxis_transform())

    ax.axhline(0, color=C_BLACK, linewidth=0.5)
    ax.set_xticks([])
    ax.set_ylabel("Composite score", fontsize=10)
    ax.set_ylim(0, max(merged[["fixed", "centered", "uniform"]].max().max() * 1.15, 1.0))
    ax.set_title("G · Three-way method comparison — Top-10 composite scores",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")
    _style_ax(ax)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_uniform_three_way_comparison.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_uniform_three_way_comparison.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def _domain_group(domains: str) -> str:
    if not isinstance(domains, str) or not domains.strip():
        return "none"
    d = domains.lower()
    if "bhlh" in d:
        return "bHLH"
    if "homeobox" in d:
        return "Homeobox"
    if "znf" in d:
        return "Znf"
    if "fork_head" in d:
        return "fork_head"
    if "t-box" in d:
        return "T-box"
    if "ets" in d:
        return "Ets"
    return "none"


def _add_panel_letter(ax, letter, x=-0.08, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top", ha="left", color=C_BLACK)


def _style_ax(ax):
    ax.xaxis.grid(True, color=C_GRID, linewidth=0.5, zorder=0)
    ax.yaxis.grid(True, color=C_GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(colors=C_BLACK)
    for spine in ax.spines.values():
        spine.set_color(C_BLACK)
    return ax


def main() -> int:
    print("== Visualize Dirichlet-uniform (alpha=1) prioritization ==")
    fig_trackA_top5()
    fig_trackB_top5()
    fig_scatter_fixed_vs_uniform()
    fig_scatter_centered_vs_uniform()
    fig_combined_dual_track()
    fig_score_shift()
    fig_three_way_comparison()
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())