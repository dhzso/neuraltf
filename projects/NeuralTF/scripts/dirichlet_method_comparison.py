#!/usr/bin/env python
"""3-way method comparison figures for Dirichlet sensitivity analysis.

Reads (from `projects/NeuralTF/results/`, gitignored):
  - dirichlet_top10_prioritized.csv          (centered)
  - dirichlet_uniform_top10.csv              (uniform)
  - dirichlet_uniform_full_rank.csv          (all 99 with both scores)

Outputs (PNGs into `projects/NeuralTF/figures/`, gitignored):
  - fig_dirichlet_score_density.png          (KDE of all 99 candidates, 3 methods)
  - fig_dirichlet_rank_correlation.png       (Spearman correlation heatmap)
  - fig_dirichlet_score_volatility.png       (per-candidate score range)
  - fig_dirichlet_method_summary.png         (4-panel summary)

Usage:
    python projects/NeuralTF/scripts/dirichlet_method_comparison.py
"""
from __future__ import annotations

from pathlib import Path
import sys

# Import shared helpers
sys.path.insert(0, str(Path(__file__).resolve().parent))
from dirichlet_visualize import (
    _style_ax, _domain_group,
    C_A, C_B, C_GRAY, C_GRID, C_BLACK, DOMAIN_COLORS,
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import spearmanr

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

C_FIXED = "#666666"
C_CENTERED = "#E69F00"
C_UNIFORM = "#009E73"


# ---------------------------------------------------------------------------
# Helper: build all-99 score table for 3 methods
# ---------------------------------------------------------------------------
def _build_all_99_scores() -> pd.DataFrame:
    """Return DataFrame with fixed, centered, uniform scores for all 99 candidates."""
    RUN = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
    rank = pd.read_csv(RUN / "rank_neural.csv")
    S = rank[STREAMS].to_numpy(dtype=float)
    mask = ~np.isnan(S)
    rank["fixed_score"] = np.where(mask, S, 0.0) @ W_DEFAULT

    full = pd.read_csv(IN_DIR / "dirichlet_uniform_full_rank.csv")
    full = full.rename(columns={"uniform_median_score": "uniform_score"})

    # Get centered score from rank_neural (it was computed by dirichlet_prioritize.py
    # but not saved in rank.csv — read from overall top-10 + full rank if needed)
    # Since dirichlet_prioritize.py doesn't save a full-rank CSV, we only have
    # centered scores for candidates that appear in dirichlet_top10_prioritized.csv
    centered_csv = IN_DIR / "dirichlet_top10_prioritized.csv"
    if centered_csv.exists():
        centered = pd.read_csv(centered_csv)
        centered_map = dict(zip(centered["gene_id_v6"], centered["dirichlet_median_score"]))
        rank["centered_score"] = rank["gene_id"].map(centered_map)
    else:
        rank["centered_score"] = np.nan

    df = rank[["gene_id", "gene_name", "fixed_score", "centered_score"]].copy()
    df = df.merge(full[["gene_id_v6", "uniform_score"]],
                  left_on="gene_id", right_on="gene_id_v6", how="left")
    df = df.drop(columns=["gene_id_v6"])
    return df


# ---------------------------------------------------------------------------
# Figure 1: Score density (KDE) for all 99 candidates, 3 methods
# ---------------------------------------------------------------------------
def fig_score_density() -> None:
    df = _build_all_99_scores()

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    x_grid = np.linspace(0, 1.0, 200)

    for col, color, label in [
        ("fixed_score", C_FIXED, "Fixed-weight (W defaults)"),
        ("centered_score", C_CENTERED, "Centered Dirichlet (k=40, all 99)"),
        ("uniform_score", C_UNIFORM, "Uniform Dirichlet (α=1, all 99)"),
    ]:
        vals = df[col].dropna().values
        if len(vals) < 2:
            continue
        # KDE
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(vals, bw_method=0.3)
        density = kde(x_grid)
        ax.plot(x_grid, density, color=color, linewidth=2.0, label=label, zorder=3)
        ax.fill_between(x_grid, density, alpha=0.12, color=color, zorder=2)

    # Mark top-10 thresholds (approx: 10th highest score)
    for col, color in [("fixed_score", C_FIXED),
                        ("centered_score", C_CENTERED),
                        ("uniform_score", C_UNIFORM)]:
        vals = df[col].dropna().sort_values(ascending=False)
        if len(vals) >= 10:
            threshold = vals.iloc[9]
            ax.axvline(threshold, color=color, linestyle=":", linewidth=0.8, alpha=0.5,
                       zorder=1)

    ax.set_xlabel("Base integrated score (0–1, before composite bonuses)", fontsize=10)
    ax.set_ylabel("Density (KDE)", fontsize=10)
    ax.set_title("Score distribution across 3 methods — all 99 candidates\n"
                 "(dotted verticals = top-10 threshold for each method)",
                 fontsize=11, fontweight="bold", pad=12)
    ax.set_xlim(0, 1.0)
    ax.legend(loc="upper right", fontsize=8, title="Method", title_fontsize=8)
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_score_density.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_score_density.png")


# ---------------------------------------------------------------------------
# Figure 2: Rank correlation matrix (Spearman) between 3 methods
# ---------------------------------------------------------------------------
def fig_rank_correlation() -> None:
    df = _build_all_99_scores()
    methods = ["fixed_score", "centered_score", "uniform_score"]
    labels = ["Fixed-weight", "Centered Dirichlet", "Uniform Dirichlet"]
    n = len(methods)

    corr = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            sub = df[[methods[i], methods[j]]].dropna()
            if len(sub) >= 3:
                result = spearmanr(sub[methods[i]], sub[methods[j]])
                # scipy.stats.spearmanr returns (rho, pvalue) for 1D inputs
                rho = result.statistic if hasattr(result, "statistic") else result[0]
                # If still a 2D matrix (2D inputs), extract diagonal
                if hasattr(rho, "__len__") and not np.isscalar(rho):
                    rho = float(np.asarray(rho).flat[0])
                corr[i, j] = float(rho)

    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    im = ax.imshow(corr, cmap="RdYlGn", vmin=0.5, vmax=1.0, zorder=2)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=10, fontweight="bold")
    ax.set_yticklabels(labels, fontsize=10, fontweight="bold")

    for i in range(n):
        for j in range(n):
            if not np.isnan(corr[i, j]):
                color = "white" if corr[i, j] < 0.85 else C_BLACK
                ax.text(j, i, f"{corr[i, j]:.3f}", ha="center", va="center",
                        fontsize=11, fontweight="bold", color=color, zorder=3)

    # Note: centered vs uniform uses subset of candidates
    n_centered = df["centered_score"].notna().sum()
    n_uniform = df["uniform_score"].notna().sum()
    n_both = (df["centered_score"].notna() & df["uniform_score"].notna()).sum()

    ax.set_title(f"Spearman rank correlation between methods\n"
                 f"(centered: {n_centered}/99, uniform: {n_uniform}/99, overlap: {n_both}/99)",
                 fontsize=11, fontweight="bold", pad=12)
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Spearman ρ", fontsize=10)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_rank_correlation.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_rank_correlation.png")


# ---------------------------------------------------------------------------
# Figure 3: Score volatility — per-candidate range across 3 methods
# ---------------------------------------------------------------------------
def fig_score_volatility() -> None:
    df = _build_all_99_scores()
    # Use only candidates that have all 3 scores
    df_complete = df.dropna(subset=["fixed_score", "centered_score", "uniform_score"])
    df_complete = df_complete.copy()
    df_complete["max_score"] = df_complete[["fixed_score", "centered_score",
                                             "uniform_score"]].max(axis=1)
    df_complete["min_score"] = df_complete[["fixed_score", "centered_score",
                                             "uniform_score"]].min(axis=1)
    df_complete["range"] = df_complete["max_score"] - df_complete["min_score"]
    df_complete = df_complete.sort_values("range", ascending=False).head(20)

    fig, ax = plt.subplots(figsize=(7.5, 6.0))
    y = np.arange(len(df_complete))
    fixed = df_complete["fixed_score"].values
    centered = df_complete["centered_score"].values
    uniform = df_complete["uniform_score"].values

    # Draw lines connecting min to max
    for i, row in df_complete.iterrows():
        x_min = row["min_score"]
        x_max = row["max_score"]
        ax.plot([x_min, x_max], [y[df_complete.index.get_loc(i)], y[df_complete.index.get_loc(i)]],
                color=C_GRAY, linewidth=1.0, zorder=2)

    # Plot the 3 scores per candidate
    ax.scatter(fixed, y, color=C_FIXED, s=50, zorder=3, label="Fixed", alpha=0.9)
    ax.scatter(centered, y, color=C_CENTERED, s=50, zorder=3, label="Centered", alpha=0.9)
    ax.scatter(uniform, y, color=C_UNIFORM, s=50, zorder=3, label="Uniform", alpha=0.9)

    ax.set_yticks(y)
    ax.set_yticklabels(df_complete["gene_name"], fontsize=9, fontweight="bold")
    ax.set_xlabel("Score range (min → max across 3 methods)", fontsize=10)
    ax.set_xlim(0, 1.0)
    ax.set_title("Score volatility across 3 methods — top 20 candidates by range\n"
                 "(candidates with all 3 scores; candidates with high range = "
                 "weight-sensitive)",
                 fontsize=11, fontweight="bold", pad=12)
    ax.legend(loc="lower right", fontsize=8)
    _style_ax(ax)

    # Add range annotation on the right
    for i, (_, row) in enumerate(df_complete.iterrows()):
        ax.text(0.98, i, f"Δ={row['range']:.3f}", transform=ax.get_yaxis_transform(),
                va="center", ha="right", fontsize=7, color=C_BLACK)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_score_volatility.png")
    plt.close(fig)
    print("  wrote fig_dirichlet_score_volatility.png")


# ---------------------------------------------------------------------------
# Figure 4: 4-panel summary (compact comparison)
# ---------------------------------------------------------------------------
def fig_method_summary() -> None:
    """Compact 4-panel summary of all 3 methods.

    Panel A: top-10 by each method (overlap shown as a Venn-like bar)
    Panel B: track-by-track winner distribution
    Panel C: score volatility distribution (histogram)
    Panel D: composite vs base score scatter (shows bonus effect)
    """
    centered = pd.read_csv(IN_DIR / "dirichlet_top10_prioritized.csv")
    uniform = pd.read_csv(IN_DIR / "dirichlet_uniform_top10.csv")
    fixed_top10 = pd.read_csv(REPO / "projects" / "NeuralTF" / "results"
                              / "top10_neural_tfs_prioritized.csv")

    fixed_ids = set(fixed_top10["gene_id_v6"])
    centered_ids = set(centered["gene_id_v6"])
    uniform_ids = set(uniform["gene_id_v6"])

    df = _build_all_99_scores()
    df_complete = df.dropna(subset=["fixed_score", "centered_score", "uniform_score"]).copy()
    df_complete["range"] = (df_complete[["fixed_score", "centered_score",
                                          "uniform_score"]].max(axis=1)
                           - df_complete[["fixed_score", "centered_score",
                                          "uniform_score"]].min(axis=1))

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    # Panel A: Venn-like bar
    ax = axes[0, 0]
    all_three = fixed_ids & centered_ids & uniform_ids
    fixed_only = fixed_ids - centered_ids - uniform_ids
    centered_only = centered_ids - fixed_ids - uniform_ids
    uniform_only = uniform_ids - fixed_ids - centered_ids
    f_c = (fixed_ids & centered_ids) - uniform_ids
    f_u = (fixed_ids & uniform_ids) - centered_ids
    c_u = (centered_ids & uniform_ids) - fixed_ids

    categories = ["All 3", "Fixed∩Center\n(−Uniform)", "Fixed∩Uniform\n(−Centered)",
                  "Center∩Uniform\n(−Fixed)", "Fixed only", "Centered only", "Uniform only"]
    sizes = [len(all_three), len(f_c), len(f_u), len(c_u),
             len(fixed_only), len(centered_only), len(uniform_only)]
    colors_a = [C_BLACK, "#888888", "#bbbbbb", "#dddddd",
                C_FIXED, C_CENTERED, C_UNIFORM]

    bars = ax.barh(categories, sizes, color=colors_a, edgecolor="white", zorder=3)
    ax.set_xlabel("Number of candidates", fontsize=10)
    ax.set_title("A · Top-10 overlap between methods", fontsize=11, fontweight="bold")
    ax.set_xlim(0, 11)
    for i, (cat, sz) in enumerate(zip(categories, sizes)):
        if sz > 0:
            ax.text(sz + 0.2, i, str(sz), va="center", fontsize=9, fontweight="bold")
    _style_ax(ax)
    ax.invert_yaxis()

    # Panel B: Track distribution in top-10 of each method
    ax = axes[0, 1]
    methods_b = ["Fixed", "Centered", "Uniform"]
    a_counts = [len(fixed_top10[fixed_top10["track"] == "A"]),
                len(centered[centered["track"] == "A"]),
                len(uniform[uniform["track"] == "A"])]
    b_counts = [len(fixed_top10[fixed_top10["track"] == "B"]),
                len(centered[centered["track"] == "B"]),
                len(uniform[uniform["track"] == "B"])]

    x = np.arange(len(methods_b))
    width = 0.35
    ax.bar(x - width/2, a_counts, width, color=C_A, edgecolor="white",
           label="Track A (RNAi-validated)", zorder=3)
    ax.bar(x + width/2, b_counts, width, color=C_B, edgecolor="white",
           label="Track B (novel)", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(methods_b, fontsize=10, fontweight="bold")
    ax.set_ylabel("Number of candidates in top-10", fontsize=10)
    ax.set_ylim(0, 6)
    ax.set_title("B · Track composition of top-10", fontsize=11, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8)
    for i, (a, b) in enumerate(zip(a_counts, b_counts)):
        ax.text(i - width/2, a + 0.1, str(a), ha="center", fontsize=9, fontweight="bold")
        ax.text(i + width/2, b + 0.1, str(b), ha="center", fontsize=9, fontweight="bold")
    _style_ax(ax)

    # Panel C: Score volatility histogram
    ax = axes[1, 0]
    ax.hist(df_complete["range"], bins=15, color=C_GRAY, edgecolor="white", zorder=3)
    median_range = df_complete["range"].median()
    ax.axvline(median_range, color=C_BLACK, linestyle="--", linewidth=1.2,
               label=f"median = {median_range:.3f}")
    ax.set_xlabel("Score range across 3 methods (max − min)", fontsize=10)
    ax.set_ylabel("Number of candidates", fontsize=10)
    ax.set_title(f"C · Volatility distribution (n={len(df_complete)})",
                 fontsize=11, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8)
    _style_ax(ax)

    # Panel D: Composite vs base scatter (shows bonus effect)
    ax = axes[1, 1]
    # Use uniform scores + composite from the uniform top-10 CSV
    for _, r in uniform.iterrows():
        gid = r["gene_id_v6"]
        if gid in df_complete["gene_id"].values:
            base = df_complete[df_complete["gene_id"] == gid]["uniform_score"].iloc[0]
            comp = r["composite_score"]
            color = C_A if r["track"] == "A" else C_B
            ax.scatter(base, comp, color=color, s=60, edgecolors="white", zorder=3)
            ax.annotate(r["gene_name"], (base, comp),
                        textcoords="offset points", xytext=(5, 3),
                        fontsize=7, color=C_BLACK)
    ax.plot([0, 1], [0, 1], "--", color=C_GRAY, linewidth=1.0, zorder=1,
            label="y = x (no bonus)")
    ax.set_xlabel("Base uniform-Dirichlet score", fontsize=10)
    ax.set_ylabel("Composite score (base + bonuses)", fontsize=10)
    ax.set_xlim(0, 0.85)
    ax.set_ylim(0, 1.0)
    ax.set_title("D · Composite bonus effect (uniform top-10)",
                 fontsize=11, fontweight="bold")
    handles = [
        Patch(facecolor=C_A, label="Track A (RNAi-validated)"),
        Patch(facecolor=C_B, label="Track B (novel)"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8)
    _style_ax(ax)

    fig.suptitle("3-way method comparison summary — Dirichlet sensitivity analysis",
                 fontsize=13, fontweight="bold", y=1.00)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_method_summary.png",
                bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_dirichlet_method_summary.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print("== 3-way method comparison figures ==")
    fig_score_density()
    fig_rank_correlation()
    fig_score_volatility()
    fig_method_summary()
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
