#!/usr/bin/env python
"""Publication-quality figures for Dirichlet-UNIFORM prioritization (Nature style).

Same 5 figures as `dirichlet_visualize.py` but for the uniform-prior analysis.
Mirrors the centered-Dirichlet visualization to allow direct comparison.

Reads (from `projects/NeuralTF/results/`, gitignored):
  - dirichlet_uniform_top10.csv         (5 Track A + 5 Track B)
  - dirichlet_uniform_full_rank.csv     (all 99 candidates with both scores)

Outputs (PNGs into `projects/NeuralTF/figures/`, gitignored):
  - fig_dirichlet_uniform_trackA_top5.png
  - fig_dirichlet_uniform_trackB_top5.png
  - fig_dirichlet_uniform_scatter.png     (fixed vs uniform)
  - fig_dirichlet_uniform_combined.png    (Track A + B with composite)
  - fig_dirichlet_uniform_score_shift.png (uniform - fixed weight)

Usage:
    python projects/NeuralTF/scripts/dirichlet_uniform_full_figures.py
"""
from __future__ import annotations

from pathlib import Path
import sys

# Import shared helpers from the centered-Dirichlet visualization
sys.path.insert(0, str(Path(__file__).resolve().parent))
from dirichlet_visualize import (
    _style_ax, _add_panel_letter, _domain_group,
    C_A, C_B, C_GRAY, C_GRID, C_BLACK, DOMAIN_COLORS,
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO = Path(__file__).resolve().parents[3]
IN_DIR  = REPO / "projects" / "NeuralTF" / "results"
OUT_DIR = REPO / "projects" / "NeuralTF" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

STREAMS = ["expression", "specificity", "reproducibility", "rnai",
           "correlation", "neural_enriched", "neural_specificity"]
W_DEFAULT = np.array([0.211, 0.105, 0.158, 0.158, 0.105, 0.158, 0.105])


# ---------------------------------------------------------------------------
# Helper: compute fixed-weight scores from rank_neural.csv
# ---------------------------------------------------------------------------
def _fixed_score_map() -> dict[str, float]:
    RUN = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
    rank = pd.read_csv(RUN / "rank_neural.csv")
    S = rank[STREAMS].to_numpy(dtype=float)
    mask = ~np.isnan(S)
    fixed_scores = np.where(mask, S, 0.0) @ W_DEFAULT
    return dict(zip(rank["gene_id"], fixed_scores))


# ---------------------------------------------------------------------------
# Figure 1: Track A top-5
# ---------------------------------------------------------------------------
def fig_uniform_trackA() -> None:
    top = pd.read_csv(IN_DIR / "dirichlet_uniform_top10.csv")
    a = top[top["track"] == "A"].copy()
    if a.empty:
        return
    a = a.sort_values("uniform_median_score", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(7.8, 3.6))
    _add_panel_letter(ax, "A")
    y = np.arange(len(a))
    dom_groups = [_domain_group(d) for d in a["interpro_domains"]]
    colors = [DOMAIN_COLORS.get(g, C_GRAY) for g in dom_groups]

    bars = ax.barh(y, a["uniform_median_score"], color=colors,
                   edgecolor="white", linewidth=0.6, height=0.6, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(a["gene_name"], fontsize=10, fontweight="bold")
    ax.set_xlabel("Uniform-Dirichlet median integrated score (α=1)", fontsize=10)
    ax.set_xlim(0, 0.85)
    ax.set_title("Track A — RNAi-validated neural TFs (uniform prior)",
                 fontsize=11, fontweight="bold", pad=12)
    _style_ax(ax)

    for i, (_, r) in enumerate(a.iterrows()):
        dom = dom_groups[i]
        comp = r["composite_score"]
        ax.text(r["uniform_median_score"] + 0.012, i,
                f'{r["uniform_median_score"]:.3f}  (composite {comp:.3f})  [{dom}]',
                va="center", fontsize=8.5, color=C_BLACK)

    present = sorted(set(dom_groups), key=lambda g: ["bHLH", "Homeobox", "Znf", "fork_head", "T-box", "Ets", "none"].index(g) if g in ["bHLH", "Homeobox", "Znf", "fork_head", "T-box", "Ets", "none"] else 99)
    handles = [Line2D([0], [0], marker="s", linestyle="", markersize=8,
                       color=DOMAIN_COLORS[g], label=g) for g in present]
    ax.legend(handles=handles, loc="lower right", frameon=True,
              framealpha=0.95, edgecolor="0.8", fontsize=8,
              title="DNA-binding domain", title_fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_uniform_trackA_top5.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_uniform_trackA_top5.png")


# ---------------------------------------------------------------------------
# Figure 2: Track B top-5
# ---------------------------------------------------------------------------
def fig_uniform_trackB() -> None:
    top = pd.read_csv(IN_DIR / "dirichlet_uniform_top10.csv")
    b = top[top["track"] == "B"].copy()
    if b.empty:
        return
    b = b.sort_values("uniform_median_score", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(7.8, 3.6))
    _add_panel_letter(ax, "B")
    y = np.arange(len(b))
    dom_groups = [_domain_group(d) for d in b["interpro_domains"]]
    colors = [DOMAIN_COLORS.get(g, C_GRAY) for g in dom_groups]

    bars = ax.barh(y, b["uniform_median_score"], color=colors,
                   edgecolor="white", linewidth=0.6, height=0.6, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(b["gene_name"], fontsize=10, fontweight="bold")
    ax.set_xlabel("Uniform-Dirichlet median integrated score (α=1)", fontsize=10)
    ax.set_xlim(0, 0.85)
    ax.set_title("Track B — Novel neural-fate candidates (uniform prior)",
                 fontsize=11, fontweight="bold", pad=12)
    _style_ax(ax)

    for i, (_, r) in enumerate(b.iterrows()):
        dom = dom_groups[i]
        comp = r["composite_score"]
        ax.text(r["uniform_median_score"] + 0.012, i,
                f'{r["uniform_median_score"]:.3f}  (composite {comp:.3f})  [{dom}]',
                va="center", fontsize=8.5, color=C_BLACK)

    present = sorted(set(dom_groups), key=lambda g: ["bHLH", "Homeobox", "Znf", "fork_head", "T-box", "Ets", "none"].index(g) if g in ["bHLH", "Homeobox", "Znf", "fork_head", "T-box", "Ets", "none"] else 99)
    handles = [Line2D([0], [0], marker="s", linestyle="", markersize=8,
                       color=DOMAIN_COLORS[g], label=g) for g in present]
    ax.legend(handles=handles, loc="lower right", frameon=True,
              framealpha=0.95, edgecolor="0.8", fontsize=8,
              title="DNA-binding domain", title_fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_uniform_trackB_top5.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_uniform_trackB_top5.png")


# ---------------------------------------------------------------------------
# Figure 3: Uniform median vs fixed-weight scatter (track-based 5A+5B)
# ---------------------------------------------------------------------------
def fig_uniform_scatter() -> None:
    top = pd.read_csv(IN_DIR / "dirichlet_uniform_top10.csv")
    fixed_map = _fixed_score_map()
    top["fixed_score"] = top["gene_id_v6"].map(fixed_map)

    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    _add_panel_letter(ax, "C")
    for trk, color, label in [("A", C_A, "Track A (RNAi-validated)"),
                               ("B", C_B, "Track B (novel)")]:
        sub = top[top["track"] == trk]
        if sub.empty:
            continue
        ax.scatter(sub["fixed_score"], sub["uniform_median_score"],
                   c=color, s=70, edgecolors="white", linewidths=0.8,
                   zorder=3, label=label, alpha=0.9)
        for _, r in sub.iterrows():
            ax.annotate(r["gene_name"],
                        (r["fixed_score"], r["uniform_median_score"]),
                        textcoords="offset points", xytext=(6, 5),
                        fontsize=8, color=C_BLACK, fontweight="medium")

    lims = [0.30, 0.82]
    ax.plot(lims, lims, "--", color=C_GRAY, linewidth=1.0, zorder=1,
            label="y = x (no change)")
    x_grid = np.linspace(lims[0], lims[1], 100)
    ax.fill_between(x_grid, x_grid - 0.03, x_grid + 0.03,
                    color=C_GRAY, alpha=0.1, zorder=0, label="±0.03 band")

    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Fixed-weight integrated score", fontsize=10)
    ax.set_ylabel("Uniform-Dirichlet median score (α=1)", fontsize=10)
    ax.set_title("Uniform-prior robustness: fixed vs uniform-Dirichlet scoring\n"
                 "(Track-based 5A + 5B)",
                 fontsize=11, fontweight="bold", pad=12)
    ax.legend(loc="upper left", fontsize=8)
    _style_ax(ax)
    ax.set_aspect('equal', adjustable='box')
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_uniform_scatter.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_uniform_scatter.png")


# ---------------------------------------------------------------------------
# Figure 4: Combined Track A + B with composite bonus overlay
# ---------------------------------------------------------------------------
def fig_uniform_combined() -> None:
    top = pd.read_csv(IN_DIR / "dirichlet_uniform_top10.csv")
    a = top[top["track"] == "A"].sort_values("uniform_median_score", ascending=False)
    b = top[top["track"] == "B"].sort_values("uniform_median_score", ascending=False)
    combined = pd.concat([a, b], ignore_index=True)

    fig, ax = plt.subplots(figsize=(9.5, 5.0))
    _add_panel_letter(ax, "D")
    x = np.arange(len(combined))
    width = 0.35
    colors = [C_A if t == "A" else C_B for t in combined["track"]]

    bars1 = ax.bar(x - width / 2, combined["uniform_median_score"], width,
                   color=colors, edgecolor="white", linewidth=0.5,
                   label="Uniform median score", zorder=3, alpha=0.95)
    bars2 = ax.bar(x + width / 2, combined["composite_score"], width,
                   color=colors, edgecolor="white", linewidth=0.5,
                   label="Composite score (with bonuses)", zorder=3,
                   alpha=0.45, hatch="///")

    ax.set_xticks(x)
    ax.set_xticklabels(combined["gene_name"], rotation=30, ha="right",
                       fontsize=9, fontweight="bold")
    ax.set_ylabel("Score", fontsize=10)
    ax.set_ylim(0, 1.10)
    ax.set_title("Uniform-Dirichlet dual-track prioritization (α=1)",
                 fontsize=12, fontweight="bold", pad=16)
    _style_ax(ax)

    ax.axvline(len(a) - 0.5, color=C_GRAY, linestyle=":", linewidth=1.0, zorder=2)

    for i, (_, r) in enumerate(combined.iterrows()):
        dom = _domain_group(r["interpro_domains"])
        ax.text(x[i] + width / 2,
                max(r["uniform_median_score"], r["composite_score"]) + 0.03,
                f"[{dom}]", ha="center", va="bottom", fontsize=7,
                color=C_BLACK, rotation=0)

    handles = [
        Line2D([0], [0], marker="s", linestyle="", markersize=9,
               color=C_A, label="Track A (RNAi-validated)"),
        Line2D([0], [0], marker="s", linestyle="", markersize=9,
               color=C_B, label="Track B (novel)"),
        Line2D([0], [0], marker="s", linestyle="", markersize=9,
               color="0.4", label="Uniform median"),
        Line2D([0], [0], marker="s", linestyle="", markersize=9,
               color="0.4", alpha=0.45, label="Composite (+ bonuses)"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8,
              title="Legend", title_fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_uniform_combined.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_uniform_combined.png")


# ---------------------------------------------------------------------------
# Figure 5: Score shift (Uniform - fixed) for track-based top-10 (5A+5B)
# ---------------------------------------------------------------------------
def fig_uniform_score_shift() -> None:
    top = pd.read_csv(IN_DIR / "dirichlet_uniform_top10.csv")
    fixed_map = _fixed_score_map()
    top["fixed_score"] = top["gene_id_v6"].map(fixed_map)
    top["diff"] = top["uniform_median_score"] - top["fixed_score"]
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
    ax.set_xlabel("Uniform-Dirichlet median − Fixed-weight score", fontsize=10)
    ax.set_title("Score shift under uniform-Dirichlet sampling\n"
                 "(5 Track A + 5 Track B)",
                 fontsize=11, fontweight="bold", pad=12)
    _style_ax(ax)

    max_abs_diff = max(abs(top["diff"].min()), abs(top["diff"].max()))
    padding = max_abs_diff * 0.3 + 0.005
    ax.set_xlim(top["diff"].min() - padding, top["diff"].max() + padding)

    for i, (_, r) in enumerate(top.iterrows()):
        val = r["diff"]
        if pd.isna(val):
            continue
        if abs(val) > 0.015:
            ax.text(val / 2, i, f"{val:+.4f}", va="center", ha="center",
                    fontsize=8, color="white", fontweight="bold")
        else:
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
    fig.savefig(OUT_DIR / "fig_dirichlet_uniform_score_shift.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_uniform_score_shift.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print("== Uniform-Dirichlet prioritization figures (Nature style) ==")
    fig_uniform_trackA()
    fig_uniform_trackB()
    fig_uniform_scatter()
    fig_uniform_combined()
    fig_uniform_score_shift()
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
