#!/usr/bin/env python
"""Publication-quality figures for Dirichlet-robust prioritization (Nature style).

Reads (from `projects/NeuralTF/results/`, gitignored):
  - dirichlet_top10_prioritized.csv   (5 Track A + 5 Track B)
  - dirichlet_overall_top10.csv        (overall top-10 by Dirichlet median)

Outputs (PNGs into `projects/NeuralTF/figures/`, gitignored):
  - fig_dirichlet_trackA_top5.png
  - fig_dirichlet_trackB_top5.png
  - fig_dirichlet_scatter.png
  - fig_dirichlet_combined.png
  - fig_dirichlet_score_shift.png

Usage:
    python projects/NeuralTF/scripts/dirichlet_visualize.py
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

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO = Path(__file__).resolve().parents[3]
IN_DIR  = REPO / "projects" / "NeuralTF" / "results"
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
C_A = "#D55E00"   # vermillion — Track A (RNAi-validated)
C_B = "#0072B2"   # blue       — Track B (novel)
C_GRAY = "#999999"
C_GRID = "#E0E0E0"
C_BLACK = "#111111"
C_GREEN = "#009E73"
C_SKY = "#56B4E9"

DOMAIN_COLORS = {
    "bHLH": "#0072B2",
    "Homeobox": "#009E73",
    "Znf": "#CC79A7",
    "fork_head": "#E69F00",
    "T-box": "#D55E00",
    "Ets": "#56B4E9",
    "none": "#BBBBBB",
}


def _domain_group(domains: str) -> str:
    d = str(domains).lower()
    if "bhlh" in d:
        return "bHLH"
    if "homeobox" in d or "homeodomain" in d:
        return "Homeobox"
    if "znf" in d or "c2h2" in d:
        return "Znf"
    if "fork_head" in d or "forkhead" in d:
        return "fork_head"
    if "t-box" in d or "t_box" in d:
        return "T-box"
    if "ets" in d:
        return "Ets"
    return "none"


def _style_ax(ax):
    ax.xaxis.grid(True, color=C_GRID, linewidth=0.5, zorder=0)
    ax.yaxis.grid(True, color=C_GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(colors=C_BLACK)
    for spine in ax.spines.values():
        spine.set_color(C_BLACK)
    return ax


def _add_panel_letter(ax, letter: str):
    ax.text(-0.18, 1.06, letter, transform=ax.transAxes,
            fontsize=12, fontweight="bold", fontstyle="italic")


# ---------------------------------------------------------------------------
# Figure 1: Track A top-5
# ---------------------------------------------------------------------------
def fig_trackA() -> None:
    top = pd.read_csv(IN_DIR / "dirichlet_top10_prioritized.csv")
    a = top[top["track"] == "A"].copy()
    if a.empty:
        return
    a = a.sort_values("dirichlet_median_score", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(7.8, 3.6))
    _add_panel_letter(ax, "A")
    y = np.arange(len(a))
    dom_groups = [_domain_group(d) for d in a["interpro_domains"]]
    colors = [DOMAIN_COLORS.get(g, C_GRAY) for g in dom_groups]

    bars = ax.barh(y, a["dirichlet_median_score"], color=colors,
                   edgecolor="white", linewidth=0.6, height=0.6, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(a["gene_name"], fontsize=10, fontweight="bold")
    ax.set_xlabel("Dirichlet median integrated score", fontsize=10)
    ax.set_xlim(0, 0.85)
    ax.set_title("Track A — RNAi-validated neural TFs",
                 fontsize=11, fontweight="bold", pad=12)
    _style_ax(ax)

    for i, (_, r) in enumerate(a.iterrows()):
        dom = dom_groups[i]
        comp = r["composite_score"]
        ax.text(r["dirichlet_median_score"] + 0.012, i,
                f'{r["dirichlet_median_score"]:.3f}  (composite {comp:.3f})  [{dom}]',
                va="center", fontsize=8.5, color=C_BLACK)

    present = sorted(set(dom_groups), key=lambda g: ["bHLH", "Homeobox", "Znf", "fork_head", "T-box", "Ets", "none"].index(g) if g in ["bHLH", "Homeobox", "Znf", "fork_head", "T-box", "Ets", "none"] else 99)
    handles = [Line2D([0], [0], marker="s", linestyle="", markersize=8,
                       color=DOMAIN_COLORS[g], label=g) for g in present]
    ax.legend(handles=handles, loc="lower right", frameon=True,
              framealpha=0.95, edgecolor="0.8", fontsize=8,
              title="DNA-binding domain", title_fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_trackA_top5.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_trackA_top5.png")


# ---------------------------------------------------------------------------
# Figure 2: Track B top-5
# ---------------------------------------------------------------------------
def fig_trackB() -> None:
    top = pd.read_csv(IN_DIR / "dirichlet_top10_prioritized.csv")
    b = top[top["track"] == "B"].copy()
    if b.empty:
        return
    b = b.sort_values("dirichlet_median_score", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(7.8, 3.6))
    _add_panel_letter(ax, "B")
    y = np.arange(len(b))
    dom_groups = [_domain_group(d) for d in b["interpro_domains"]]
    colors = [DOMAIN_COLORS.get(g, C_GRAY) for g in dom_groups]

    bars = ax.barh(y, b["dirichlet_median_score"], color=colors,
                   edgecolor="white", linewidth=0.6, height=0.6, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(b["gene_name"], fontsize=10, fontweight="bold")
    ax.set_xlabel("Dirichlet median integrated score", fontsize=10)
    ax.set_xlim(0, 0.85)
    ax.set_title("Track B — Novel neural-fate candidates",
                 fontsize=11, fontweight="bold", pad=12)
    _style_ax(ax)

    for i, (_, r) in enumerate(b.iterrows()):
        dom = dom_groups[i]
        comp = r["composite_score"]
        ax.text(r["dirichlet_median_score"] + 0.012, i,
                f'{r["dirichlet_median_score"]:.3f}  (composite {comp:.3f})  [{dom}]',
                va="center", fontsize=8.5, color=C_BLACK)

    present = sorted(set(dom_groups), key=lambda g: ["bHLH", "Homeobox", "Znf", "fork_head", "T-box", "Ets", "none"].index(g) if g in ["bHLH", "Homeobox", "Znf", "fork_head", "T-box", "Ets", "none"] else 99)
    handles = [Line2D([0], [0], marker="s", linestyle="", markersize=8,
                       color=DOMAIN_COLORS[g], label=g) for g in present]
    ax.legend(handles=handles, loc="lower right", frameon=True,
              framealpha=0.95, edgecolor="0.8", fontsize=8,
              title="DNA-binding domain", title_fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_trackB_top5.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_trackB_top5.png")


# ---------------------------------------------------------------------------
# Figure 3: Dirichlet median vs fixed-weight scatter (track-based 5A+5B)
# ---------------------------------------------------------------------------
def fig_scatter() -> None:
    # Use track-based top-10 (5A+5B) from overall CSV
    overall = pd.read_csv(IN_DIR / "dirichlet_overall_top10.csv")
    # Compute fixed-weight scores for all 10 from source data
    REPO = Path(__file__).resolve().parents[3]
    RUN = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
    STREAMS = ["expression", "specificity", "reproducibility", "rnai",
               "correlation", "neural_enriched", "neural_specificity"]
    W_DEFAULT = np.array([0.211, 0.105, 0.158, 0.158, 0.105, 0.158, 0.105])
    
    rank = pd.read_csv(RUN / "rank_neural.csv")
    S = rank[STREAMS].to_numpy(dtype=float)
    mask = ~np.isnan(S)
    fixed_scores = np.where(mask, S, 0.0) @ W_DEFAULT
    fixed_map = dict(zip(rank["gene_id"], fixed_scores))
    
    overall["fixed_weight_score"] = overall["gene_id_v6"].map(fixed_map)
    plot_data = overall  # 5A + 5B

    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    _add_panel_letter(ax, "C")
    for trk, color, label in [("A", C_A, "Track A (RNAi-validated)"),
                               ("B", C_B, "Track B (novel)")]:
        sub = plot_data[plot_data["track"] == trk]
        if sub.empty:
            continue
        ax.scatter(sub["fixed_weight_score"], sub["dirichlet_median_score"],
                   c=color, s=70, edgecolors="white", linewidths=0.8,
                   zorder=3, label=label, alpha=0.9)
        for _, r in sub.iterrows():
            ax.annotate(r["gene_name"],
                        (r["fixed_weight_score"], r["dirichlet_median_score"]),
                        textcoords="offset points", xytext=(6, 5),
                        fontsize=8, color=C_BLACK, fontweight="medium")

    # Identity line
    lims = [0.30, 0.82]
    ax.plot(lims, lims, "--", color=C_GRAY, linewidth=1.0, zorder=1,
            label="y = x (no change)")

    # Confidence band
    x_grid = np.linspace(lims[0], lims[1], 100)
    ax.fill_between(x_grid, x_grid - 0.03, x_grid + 0.03,
                    color=C_GRAY, alpha=0.1, zorder=0, label="\u00B10.03 band")

    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Fixed-weight integrated score", fontsize=10)
    ax.set_ylabel("Dirichlet median score", fontsize=10)
    ax.set_title("Robustness: Dirichlet-robust vs fixed-weight scoring\n(Track-based 5A + 5B)",
                 fontsize=11, fontweight="bold", pad=12)
    ax.legend(loc="upper left", fontsize=8)
    _style_ax(ax)
    ax.set_aspect('equal', adjustable='box')
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_scatter.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_scatter.png")


# ---------------------------------------------------------------------------
# Figure 4: Combined Track A + B (main figure)
# ---------------------------------------------------------------------------
def fig_combined() -> None:
    top = pd.read_csv(IN_DIR / "dirichlet_top10_prioritized.csv")
    a = top[top["track"] == "A"].sort_values("dirichlet_median_score", ascending=False)
    b = top[top["track"] == "B"].sort_values("dirichlet_median_score", ascending=False)
    combined = pd.concat([a, b], ignore_index=True)

    fig, ax = plt.subplots(figsize=(9.5, 5.0))
    _add_panel_letter(ax, "D")
    x = np.arange(len(combined))
    width = 0.35

    colors = [C_A if t == "A" else C_B for t in combined["track"]]

    bars1 = ax.bar(x - width / 2, combined["dirichlet_median_score"], width,
                   color=colors, edgecolor="white", linewidth=0.5,
                   label="Dirichlet median score", zorder=3, alpha=0.95)
    bars2 = ax.bar(x + width / 2, combined["composite_score"], width,
                   color=colors, edgecolor="white", linewidth=0.5,
                   label="Composite score (with bonuses)", zorder=3,
                   alpha=0.45, hatch="///")

    ax.set_xticks(x)
    ax.set_xticklabels(combined["gene_name"], rotation=30, ha="right",
                       fontsize=9, fontweight="bold")
    ax.set_ylabel("Score", fontsize=10)
    ax.set_ylim(0, 1.10)
    ax.set_title("Dirichlet-robust dual-track prioritization",
                 fontsize=12, fontweight="bold", pad=16)
    _style_ax(ax)

    # Track separator
    ax.axvline(len(a) - 0.5, color=C_GRAY, linestyle=":", linewidth=1.0, zorder=2)

    # Annotate domain on top of bars
    for i, (_, r) in enumerate(combined.iterrows()):
        dom = _domain_group(r["interpro_domains"])
        ax.text(x[i] + width / 2, max(r["dirichlet_median_score"], r["composite_score"]) + 0.03,
                f"[{dom}]", ha="center", va="bottom", fontsize=7,
                color=C_BLACK, rotation=0)

    # Custom legend
    handles = [
        Line2D([0], [0], marker="s", linestyle="", markersize=9,
               color=C_A, label="Track A (RNAi-validated)"),
        Line2D([0], [0], marker="s", linestyle="", markersize=9,
               color=C_B, label="Track B (novel)"),
        Line2D([0], [0], marker="s", linestyle="", markersize=9,
               color="0.4", label="Dirichlet median"),
        Line2D([0], [0], marker="s", linestyle="", markersize=9,
               color="0.4", alpha=0.45, label="Composite (+ bonuses)"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8,
              title="Legend", title_fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_combined.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_combined.png")


# ---------------------------------------------------------------------------
# Figure 5: Score difference (Dirichlet - fixed) for track-based top-10 (5A+5B)
# ---------------------------------------------------------------------------
def fig_score_shift() -> None:
    top = pd.read_csv(IN_DIR / "dirichlet_top10_prioritized.csv")
    
    # Compute fixed-weight scores for ALL 99 candidates from rank_neural.csv
    REPO = Path(__file__).resolve().parents[3]
    RUN = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
    STREAMS = ["expression", "specificity", "reproducibility", "rnai",
               "correlation", "neural_enriched", "neural_specificity"]
    W_DEFAULT = np.array([0.211, 0.105, 0.158, 0.158, 0.105, 0.158, 0.105])
    
    rank = pd.read_csv(RUN / "rank_neural.csv")
    S = rank[STREAMS].to_numpy(dtype=float)
    mask = ~np.isnan(S)
    fixed_scores = np.where(mask, S, 0.0) @ W_DEFAULT
    fixed_map = dict(zip(rank["gene_id"], fixed_scores))
    
    top["fixed_score"] = top["gene_id_v6"].map(fixed_map)
    top["diff"] = top["dirichlet_median_score"] - top["fixed_score"]
    top = top.sort_values("diff", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    _add_panel_letter(ax, "E")
    y = np.arange(len(top))
    colors = [C_A if t == "A" else C_B for t in top["track"]]

    bars = ax.barh(y, top["diff"], color=colors,
                   edgecolor="white", linewidth=0.5, height=0.55, zorder=3)
    ax.axvline(0, color=C_BLACK, linewidth=0.8, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(top["gene_name"], fontsize=9, fontweight="bold")
    ax.set_xlabel("Dirichlet median \u2212 Fixed-weight score", fontsize=10)
    ax.set_title("Score shift under Dirichlet weight sampling\n(5 Track A + 5 Track B)",
                 fontsize=11, fontweight="bold", pad=12)
    _style_ax(ax)

    # Set xlim to accommodate text
    max_abs_diff = max(abs(top["diff"].min()), abs(top["diff"].max()))
    padding = max_abs_diff * 0.3 + 0.005
    ax.set_xlim(top["diff"].min() - padding, top["diff"].max() + padding)

    for i, (_, r) in enumerate(top.iterrows()):
        val = r["diff"]
        if pd.isna(val):
            continue
        # Place text inside bar if large enough, else outside
        if abs(val) > 0.015:
            # Inside bar
            ax.text(val / 2, i, f"{val:+.4f}", va="center", ha="center",
                    fontsize=8, color="white", fontweight="bold")
        else:
            # Outside bar
            offset = 0.002 if val >= 0 else -0.002
            ax.text(val + offset, i, f"{val:+.4f}", va="center",
                    ha="left" if val >= 0 else "right",
                    fontsize=8, color=C_BLACK, fontweight="medium")

    handles = [
        Line2D([0], [0], marker="s", linestyle="", markersize=9,
               color=C_A, label="Track A (RNAi-validated)"),
        Line2D([0], [0], marker="s", linestyle="", markersize=9,
               color=C_B, label="Track B (novel)"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_score_shift.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_score_shift.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print("== Dirichlet prioritization figures (Nature publication quality) ==")
    fig_trackA()
    fig_trackB()
    fig_scatter()
    fig_combined()
    fig_score_shift()
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())