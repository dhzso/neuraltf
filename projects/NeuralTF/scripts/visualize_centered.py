#!/usr/bin/env python
"""Visualize Dirichlet-centered (k=40) prioritization results.

Generates 5 publication-quality figures from Dirichlet-centered outputs:

  1. trackA_top5          - Track A top-5 (RNAi-validated) with domain colors
  2. trackB_top5          - Track B top-5 (novel) with domain colors
  3. scatter_fixed_vs_dirichlet - fixed-weight vs Dirichlet median
  4. combined_dual_track  - both tracks with composite overlay
  5. score_shift          - Dirichlet composite - fixed composite per candidate

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

REPO = Path(__file__).resolve().parents[3]
IN_DIR = REPO / "projects" / "NeuralTF" / "results"
OUT_DIR = REPO / "projects" / "NeuralTF" / "figures"
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
    csv_path = IN_DIR / "dirichlet_top10_prioritized.csv"
    if not csv_path.exists():
        return
    top = pd.read_csv(csv_path)
    a = top[top["track"] == "A"].copy()
    if a.empty:
        return
    a = a.sort_values("dirichlet_median_score", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(7.8, 3.6))
    _add_panel_letter(ax, "A")
    y = np.arange(len(a))
    dom_groups = [_domain_group(d) for d in a["interpro_domains"]]
    colors = [DOMAIN_COLORS.get(g, C_GRAY) for g in dom_groups]

    bars = ax.barh(y, a["composite_score"], color=colors, edgecolor="white",
                   height=0.6, linewidth=0.5, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(a["gene_name"], fontsize=9.5)
    ax.set_xlim(0, max(1.05, a["composite_score"].max() * 1.15))
    ax.set_xlabel("Composite score (Dirichlet median + bonuses)", fontsize=10)
    ax.set_title("A · Track A — RNAi-validated (Dirichlet-centered)",
                 fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    _style_ax(ax)

    # Legend
    handles = [Patch(facecolor=DOMAIN_COLORS[d], edgecolor="white",
                     label=d) for d in sorted(set(dom_groups)) if d != "none"]
    ax.legend(handles=handles, loc="lower right", fontsize=8, title="Domain group")

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

    fig, ax = plt.subplots(figsize=(7.8, 3.6))
    _add_panel_letter(ax, "B")
    y = np.arange(len(b))
    dom_groups = [_domain_group(d) for d in b["interpro_domains"]]
    colors = [DOMAIN_COLORS.get(g, C_GRAY) for g in dom_groups]

    bars = ax.barh(y, b["composite_score"], color=colors, edgecolor="white",
                   height=0.6, linewidth=0.5, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(b["gene_name"], fontsize=9.5)
    ax.set_xlim(0, max(1.05, b["composite_score"].max() * 1.15))
    ax.set_xlabel("Composite score (Dirichlet median + bonuses)", fontsize=10)
    ax.set_title("B · Track B — Novel candidates (Dirichlet-centered)",
                 fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    _style_ax(ax)

    handles = [Patch(facecolor=DOMAIN_COLORS[d], edgecolor="white",
                     label=d) for d in sorted(set(dom_groups)) if d != "none"]
    ax.legend(handles=handles, loc="lower right", fontsize=8, title="Domain group")

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

    # Merge on gene_id_v6
    merged = pd.merge(fixed[["gene_id_v6", "composite_score"]],
                      centered[["gene_id_v6", "composite_score", "dirichlet_median_score", "track"]],
                      on="gene_id_v6", suffixes=("_fixed", "_dirichlet"))

    fig, ax = plt.subplots(figsize=(6, 5))
    _add_panel_letter(ax, "C")

    for _, row in merged.iterrows():
        color = C_A if row["track"] == "A" else C_B
        marker = "o" if row["track"] == "A" else "^"
        ax.scatter(row["composite_score_fixed"], row["composite_score_dirichlet"],
                   color=color, s=80, marker=marker, edgecolors="white",
                   linewidth=0.8, zorder=3)
        ax.annotate(row["gene_id_v6"].replace("dd_Smed_v6_", "dd"),
                    (row["composite_score_fixed"], row["composite_score_dirichlet"]),
                    textcoords="offset points", xytext=(5, 3),
                    fontsize=7, color=C_BLACK)

    lims = [0, 1.0]
    ax.plot(lims, lims, "--", color=C_GRAY, linewidth=1.0, zorder=1,
            label="y = x (no shift)")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Fixed-weight composite score", fontsize=10)
    ax.set_ylabel("Dirichlet-centered composite score", fontsize=10)
    ax.set_title("C · Fixed-weight vs Dirichlet-centered composite",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
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

    fig, ax = plt.subplots(figsize=(8, 5))
    _add_panel_letter(ax, "D")

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
    ax.set_ylabel("Composite score (Dirichlet median + bonuses)", fontsize=10)
    ax.set_ylim(0, max(top["composite_score"].max() * 1.2, 1.0))
    ax.set_title("D · Dirichlet-centered dual track — Top-10",
                 fontsize=11, fontweight="bold")

    # Legends
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

    fig, ax = plt.subplots(figsize=(8, 5))
    _add_panel_letter(ax, "E")

    colors = [C_A if t == "A" else C_B for t in merged["track"]]
    x = np.arange(len(merged))
    ax.bar(x, merged["shift"], color=colors, edgecolor="white",
           width=0.6, linewidth=0.5, zorder=3)
    ax.axhline(0, color=C_BLACK, linewidth=0.8, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(merged["gene_id_v6"].str.replace("dd_Smed_v6_", "dd"),
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Dirichlet composite − Fixed composite", fontsize=10)
    ax.set_title("E · Score shift (Dirichlet-centered − fixed)",
                 fontsize=11, fontweight="bold")
    ax.legend(handles=[
        Patch(facecolor=C_A, label="Track A (RNAi)"),
        Patch(facecolor=C_B, label="Track B (novel)")
    ], loc="upper right", fontsize=8)
    _style_ax(ax)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_centered_score_shift.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_centered_score_shift.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
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