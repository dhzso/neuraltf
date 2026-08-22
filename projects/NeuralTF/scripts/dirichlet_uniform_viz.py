#!/usr/bin/env python
"""3-way comparison figure: fixed-weight vs centered-Dirichlet vs uniform-Dirichlet.

Reads (from `projects/NeuralTF/results/`, gitignored):
  - dirichlet_top10_prioritized.csv          (centered)
  - dirichlet_uniform_top10.csv              (uniform)
  - dirichlet_uniform_full_rank.csv          (all 99 with both scores)

Outputs (PNGs into `projects/NeuralTF/figures/`, gitignored):
  - fig_dirichlet_uniform_vs_centered.png    (scatter: centered vs uniform)
  - fig_dirichlet_3way_comparison.png        (grouped bars: all 3 methods)

Usage:
    python projects/NeuralTF/scripts/dirichlet_uniform_viz.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parents[3]
IN_DIR  = REPO / "projects" / "NeuralTF" / "results"
OUT_DIR = REPO / "projects" / "NeuralTF" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.12,
    "legend.frameon": True,
    "legend.framealpha": 0.92,
    "legend.edgecolor": "0.8",
})

C_A = "#D55E00"   # Track A
C_B = "#0072B2"   # Track B
C_FIXED = "#666666"
C_CENTERED = "#E69F00"
C_UNIFORM = "#009E73"
C_GRAY = "#999999"
C_GRID = "#E0E0E0"
C_BLACK = "#111111"


def _style_ax(ax):
    ax.xaxis.grid(True, color=C_GRID, linewidth=0.5, zorder=0)
    ax.yaxis.grid(True, color=C_GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(colors=C_BLACK)
    for spine in ax.spines.values():
        spine.set_color(C_BLACK)
    return ax


# ---------------------------------------------------------------------------
# Figure 1: Scatter — centered vs uniform Dirichlet
# ---------------------------------------------------------------------------
def fig_scatter_centered_vs_uniform() -> None:
    centered_csv = IN_DIR / "dirichlet_top10_prioritized.csv"
    uniform_csv = IN_DIR / "dirichlet_uniform_top10.csv"

    if not centered_csv.exists() or not uniform_csv.exists():
        print("  [skip] missing CSV")
        return

    centered = pd.read_csv(centered_csv)
    uniform = pd.read_csv(uniform_csv)

    # Merge on gene_id
    merged = centered.merge(
        uniform[["gene_id_v6", "uniform_median_score", "track"]].rename(
            columns={"track": "track_uniform"}),
        left_on="gene_id_v6", right_on="gene_id_v6", how="inner"
    )

    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    for trk, color, label in [("A", C_A, "Track A (RNAi-validated)"),
                               ("B", C_B, "Track B (novel)")]:
        sub = merged[merged["track"] == trk]
        if sub.empty:
            continue
        ax.scatter(sub["dirichlet_median_score"], sub["uniform_median_score"],
                   c=color, s=70, edgecolors="white", linewidths=0.8,
                   zorder=3, label=label, alpha=0.9)
        for _, r in sub.iterrows():
            ax.annotate(r["gene_name"],
                        (r["dirichlet_median_score"], r["uniform_median_score"]),
                        textcoords="offset points", xytext=(6, 5),
                        fontsize=8, color=C_BLACK, fontweight="medium")

    # Identity line
    lims = [0.45, 0.80]
    ax.plot(lims, lims, "--", color=C_GRAY, linewidth=1.0, zorder=1,
            label="y = x (no change)")

    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Centered-Dirichlet median (k=40, centered on defaults)",
                  fontsize=10)
    ax.set_ylabel("Uniform-Dirichlet median (α=1, no prior)", fontsize=10)
    ax.set_title("Centered vs Uniform Dirichlet median\n"
                 "(track-based 5A + 5B)",
                 fontsize=11, fontweight="bold", pad=12)
    ax.legend(loc="upper left", fontsize=8)
    _style_ax(ax)
    ax.set_aspect('equal', adjustable='box')
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_uniform_vs_centered.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_uniform_vs_centered.png")


# ---------------------------------------------------------------------------
# Figure 2: 3-way grouped bar chart
# ---------------------------------------------------------------------------
def fig_3way_comparison() -> None:
    centered_csv = IN_DIR / "dirichlet_top10_prioritized.csv"
    uniform_csv = IN_DIR / "dirichlet_uniform_top10.csv"

    if not centered_csv.exists() or not uniform_csv.exists():
        print("  [skip] missing CSV")
        return

    centered = pd.read_csv(centered_csv)
    uniform = pd.read_csv(uniform_csv)

    # Need fixed-weight scores — compute from rank_neural.csv
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

    # UNION of both track-based lists (so all candidates appear)
    all_ids = set(centered["gene_id_v6"]) | set(uniform["gene_id_v6"])
    centered_map = dict(zip(centered["gene_id_v6"], centered["dirichlet_median_score"]))
    uniform_map = dict(zip(uniform["gene_id_v6"], uniform["uniform_median_score"]))

    # Track label: take from whichever list the candidate is in
    track_map = {}
    for _, r in centered.iterrows():
        track_map[r["gene_id_v6"]] = r["track"]
    for _, r in uniform.iterrows():
        if r["gene_id_v6"] not in track_map:
            track_map[r["gene_id_v6"]] = r["track"]

    rows = []
    for gid in sorted(all_ids):
        rows.append({
            "gene_id_v6": gid,
            "gene_name": (centered[centered["gene_id_v6"] == gid].iloc[0]["gene_name"]
                         if gid in centered["gene_id_v6"].values
                         else uniform[uniform["gene_id_v6"] == gid].iloc[0]["gene_name"]),
            "track": track_map[gid],
            "fixed": fixed_map.get(gid, 0),
            "centered": centered_map.get(gid, np.nan),
            "uniform": uniform_map.get(gid, np.nan),
        })
    df = pd.DataFrame(rows)

    # Replace NaN with actual score where available — candidates not in a method's
    # top-10 still have a score from the full-rank file, so look it up.
    # Otherwise plot them at 0 with a visible "below top-10" marker.
    full_uniform_csv = IN_DIR / "dirichlet_uniform_full_rank.csv"
    if full_uniform_csv.exists():
        full_unif = pd.read_csv(full_uniform_csv)
        unif_full_map = dict(zip(full_unif["gene_id_v6"], full_unif["uniform_median_score"]))
        # Replace NaN uniform with full-rank uniform
        df["uniform"] = df.apply(
            lambda r: r["uniform"] if pd.notna(r["uniform"]) else unif_full_map.get(r["gene_id_v6"], 0.0),
            axis=1)
    # For centered: only present if in centered top-10; use fixed-score as proxy (or 0)
    df["centered"] = df["centered"].fillna(0)
    df["uniform"] = df["uniform"].fillna(0)

    # Mark which methods each candidate was a top-10 in
    df["in_centered"] = df["gene_id_v6"].isin(centered["gene_id_v6"]).astype(int)
    df["in_uniform"] = df["gene_id_v6"].isin(uniform["gene_id_v6"]).astype(int)

    # Plot: Track A on left, Track B on right, sort by fixed within track
    a = df[df["track"] == "A"].sort_values("fixed", ascending=False).reset_index(drop=True)
    b = df[df["track"] == "B"].sort_values("fixed", ascending=False).reset_index(drop=True)
    combined = pd.concat([a, b], ignore_index=True)

    fig, ax = plt.subplots(figsize=(12, 5.4))
    x = np.arange(len(combined))
    width = 0.27

    # For candidates not in a method's top-10, use a lighter shade
    def _lighten(hex_color, factor):
        h = hex_color.lstrip("#")
        rgb = np.array([int(h[i:i+2], 16) for i in (0, 2, 4)])
        white = np.array([255, 255, 255])
        mixed = (rgb * (1 - factor) + white * factor).astype(int)
        return f"#{mixed[0]:02x}{mixed[1]:02x}{mixed[2]:02x}"

    # Fixed-weight is always "in" (it's the baseline); centered/uniform use light/dark
    in_fixed = np.ones(len(combined), dtype=bool)
    in_centered = combined["in_centered"].values.astype(bool)
    in_uniform = combined["in_uniform"].values.astype(bool)

    colors_fixed = [C_FIXED if t10 else _lighten(C_FIXED, 0.55) for t10 in in_fixed]
    colors_cent = [C_CENTERED if t10 else _lighten(C_CENTERED, 0.55) for t10 in in_centered]
    colors_unif = [C_UNIFORM if t10 else _lighten(C_UNIFORM, 0.55) for t10 in in_uniform]

    bars1 = ax.bar(x - width, combined["fixed"], width,
                   color=colors_fixed, edgecolor="white", linewidth=0.5,
                   zorder=3, alpha=0.9)
    bars2 = ax.bar(x, combined["centered"], width,
                   color=colors_cent, edgecolor="white", linewidth=0.5,
                   zorder=3, alpha=0.9)
    bars3 = ax.bar(x + width, combined["uniform"], width,
                   color=colors_unif, edgecolor="white", linewidth=0.5,
                   zorder=3, alpha=0.9)

    # Add annotations for the candidates that are NOT in each method's top-10
    for i, row in combined.iterrows():
        if row["in_centered"] == 0:
            ax.text(i, 0.02, "—", ha="center", va="bottom",
                    fontsize=6, color=C_GRAY, zorder=4)
        if row["in_uniform"] == 0:
            ax.text(i + width, 0.02, "—", ha="center", va="bottom",
                    fontsize=6, color=C_GRAY, zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels(combined["gene_name"], rotation=30, ha="right",
                       fontsize=8, fontweight="bold")
    ax.set_ylabel("Base integrated score (0–1, before composite bonuses)", fontsize=10)
    ax.set_ylim(0, 1.0)
    ax.set_title("3-way score comparison: Fixed / Centered Dirichlet / Uniform Dirichlet\n"
                 "(base scores; composite bonuses applied separately for final ranking)",
                 fontsize=11, fontweight="bold", pad=14)
    _style_ax(ax)

    # Track separator
    ax.axvline(len(a) - 0.5, color=C_GRAY, linestyle=":", linewidth=1.0, zorder=2)

    # Track labels INSIDE the plot in the empty 0.8–1.0 space
    mid_a = len(a) / 2 - 0.5
    mid_b = len(a) + len(b) / 2 - 0.5
    ax.text(mid_a, 0.90, "Track A", ha="center", va="center",
            fontsize=11, fontweight="bold", color=C_A,
            bbox=dict(facecolor="white", edgecolor=C_A, boxstyle="round,pad=0.3",
                      linewidth=0.8, alpha=0.9), zorder=5)
    ax.text(mid_b, 0.90, "Track B", ha="center", va="center",
            fontsize=11, fontweight="bold", color=C_B,
            bbox=dict(facecolor="white", edgecolor=C_B, boxstyle="round,pad=0.3",
                      linewidth=0.8, alpha=0.9), zorder=5)

    # Legend (with light/dark distinction)
    handles = [
        Patch(facecolor=C_FIXED, label="Fixed-weight (top-10)"),
        Patch(facecolor=_lighten(C_FIXED, 0.55), label="Fixed-weight (below top-10)"),
        Patch(facecolor=C_CENTERED, label="Centered Dirichlet (top-10)"),
        Patch(facecolor=_lighten(C_CENTERED, 0.55), label="Centered Dirichlet (below)"),
        Patch(facecolor=C_UNIFORM, label="Uniform Dirichlet (top-10)"),
        Patch(facecolor=_lighten(C_UNIFORM, 0.55), label="Uniform Dirichlet (below)"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=7,
              title="Method", title_fontsize=8, framealpha=0.95)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_3way_comparison.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_3way_comparison.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print("== 3-way Dirichlet comparison figures ==")
    fig_scatter_centered_vs_uniform()
    fig_3way_comparison()
    print("Done.")
    return 0


# Patch: import Patch at module level
from matplotlib.patches import Patch


if __name__ == "__main__":
    raise SystemExit(main())
