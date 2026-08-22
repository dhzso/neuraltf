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
# Figure:99-neural vs 249-wide overlap (rank shift + score comparison)
# ---------------------------------------------------------------------------
def fig_99vs249() -> None:
    """Compare99-neural vs 249-wide uniform Dirichlet top-10.

    Panel A: Bar chart showing rank shift between the two scopes.
    Panel B: Score comparison (99 scope vs 249 scope scores).
    """
    # Load99-neural full rank
    unif99 = pd.read_csv(IN_DIR / "dirichlet_uniform_full_rank.csv")
    # Load249-wide full rank
    all249 = pd.read_csv(IN_DIR / "dirichlet_uniform_all249_full_rank.csv")
    # Load99-neural top-10
    unif99_top10 = pd.read_csv(IN_DIR / "dirichlet_uniform_top10.csv")
    unif99_top10_ids = set(unif99_top10["gene_id_v6"])
    # Load249-wide top-10
    all249_top10 = pd.read_csv(IN_DIR / "dirichlet_uniform_all249_top10.csv")
    all249_top10_ids = set(all249_top10["gene_id_v6"])
    # Load full rank
    RUN = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
    rank_neural = pd.read_csv(RUN / "rank_neural.csv")
    rank_all = pd.read_csv(RUN / "rank.csv")

    # Merge scores on gene_id
    scores99 = dict(zip(rank_neural["gene_id"], unif99["uniform_median_score"]))
    scores249 = dict(zip(rank_all["gene_id"], all249["uniform_median_score"]))

    # Full-rank99 list: all 99 candidates
    all99_ids = set(rank_neural["gene_id"])
    all249_ids_set = set(rank_all["gene_id"])

    # Get ranks (1-indexed, higher score = rank 1)
    rank99 = rank_neural.sort_values("integrated_score", ascending=False).reset_index(drop=True)
    rank99["rank"] = range(1, len(rank99) + 1)
    rank99_map = dict(zip(rank99["gene_id"], rank99["rank"]))

    rank249 = rank_all.sort_values("integrated_score", ascending=False).reset_index(drop=True)
    rank249["rank"] = range(1, len(rank249) + 1)
    rank249_map = dict(zip(rank249["gene_id"], rank249["rank"]))

    # Candidates in both sets (the99-neural set is a subset of 249)
    overlap_ids = all99_ids & all249_ids_set
    only_in_99 = all99_ids - all249_ids_set  # should be empty (99 is subset of 249)

    # For the overlap candidates, compute rank shift
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12, 5), dpi=150)
    fig.subplots_adjust(wspace=0.35)

    # --- Panel A: Rank shift across all99 candidates -----------------------
    rank_shift = []
    for gid in overlap_ids:
        r99 = rank99_map.get(gid, np.nan)
        r249 = rank249_map.get(gid, np.nan)
        if not np.isnan(r99) and not np.isnan(r249):
            rank_shift.append({"gene_id": gid, "rank_99": r99, "rank_249": r249,
                               "shift": r249 - r99})
    rank_shift_df = pd.DataFrame(rank_shift)

    # Color by direction of shift
    colors = []
    for s in rank_shift_df["shift"]:
        if s > 0:
            colors.append(C_B)   # rank went up (worse) in 249 = blue
        elif s < 0:
            colors.append(C_A)   # rank went down (better) in 249 = orange
        else:
            colors.append(C_GRAY) # no change

    y_pos = range(len(rank_shift_df))
    ax_a.barh(y_pos, rank_shift_df["shift"], color=colors, height=0.7,
              edgecolor="white", linewidth=0.5)
    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels(
        [rank_shift_df.iloc[i]["gene_id"].replace("dd_Smed_v6_", "dd")
         for i in range(len(rank_shift_df))],
        fontsize=7)
    ax_a.axvline(x=0, color=C_BLACK, linewidth=0.8, linestyle="--")
    ax_a.set_xlabel("Rank shift (249-wide rank − 99-neural rank)", fontsize=9)
    ax_a.set_title("A · Rank shift (99 candidates)",
                   fontsize=11, fontweight="bold")
    ax_a.invert_yaxis()
    _style_ax(ax_a)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=C_B, label="Rank higher (worse) in 249"),
        Patch(facecolor=C_A, label="Rank lower (better) in 249"),
        Patch(facecolor=C_GRAY, label="No change"),
    ]
    ax_a.legend(handles=legend_elements, loc="lower right", fontsize=7)

    # --- Panel B: Score comparison (99 scope vs 249 scope) ------------------
    # For all99 candidates, plot99_score vs 249_score
    x_scores, y_scores, labels = [], [], []
    for gid in overlap_ids:
        s99 = scores99.get(gid, np.nan)
        s249 = scores249.get(gid, np.nan)
        if not np.isnan(s99) and not np.isnan(s249):
            x_scores.append(s99)
            y_scores.append(s249)
            labels.append(gid.replace("dd_Smed_v6_", "dd"))
    x_scores = np.array(x_scores)
    y_scores = np.array(y_scores)

    # Mark top-10 candidates
    top10_99 = {g.replace("dd_Smed_v6_", "dd") for g in unif99_top10_ids}
    top10_249 = {g.replace("dd_Smed_v6_", "dd") for g in all249_top10_ids}

    for i, lbl in enumerate(labels):
        is_top10_99 = lbl in top10_99
        is_top10_249 = lbl in top10_249
        if is_top10_99 and is_top10_249:
            color = "#009E73"  # both = green
            marker = "o"
            zord = 4
        elif is_top10_99:
            color = C_CENTERED  # 99 only = orange
            marker = "^"
            zord = 3
        elif is_top10_249:
            color = C_A  # 249 only = blue
            marker = "s"
            zord = 3
        else:
            color = C_GRAY  # neither = gray
            marker = "x"
            zord = 2
        ax_b.scatter(x_scores[i], y_scores[i], color=color, s=50,
                     marker=marker, edgecolors="white", zorder=zord)
        if is_top10_99 or is_top10_249:
            ax_b.annotate(lbl, (x_scores[i], y_scores[i]),
                          textcoords="offset points", xytext=(5, 3),
                          fontsize=6, color=C_BLACK)

    # Diagonal (y=x line)
    lims = [min(x_scores.min(), y_scores.min()) - 0.02,
            max(x_scores.max(), y_scores.max()) + 0.02]
    ax_b.plot(lims, lims, "--", color=C_GRAY, linewidth=1.0, zorder=1,
              label="y = x (no change)")
    ax_b.set_xlim(lims)
    ax_b.set_ylim(lims)
    ax_b.set_xlabel("Uniform median score (99-neural scope)", fontsize=9)
    ax_b.set_ylabel("Uniform median score (249-wide scope)", fontsize=9)
    ax_b.set_title("B · Score comparison",
                   fontsize=11, fontweight="bold")
    ax_b.legend(fontsize=7)

    # Legend for markers
    legend_elements2 = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#009E73",
               markersize=8, label="Top-10 in both scopes"),
        Line2D([0], [0], marker="^", color="w", markerfacecolor=C_CENTERED,
               markersize=8, label="Top-10 only in 99-neural"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor=C_A,
               markersize=8, label="Top-10 only in 249-wide"),
        Line2D([0], [0], marker="x", color=C_GRAY, markersize=8,
               label="Not in either top-10"),
    ]
    ax_b.legend(handles=legend_elements2, loc="lower right", fontsize=7)
    _style_ax(ax_b)

    fig.suptitle("99-neural vs 249-wide Dirichlet-uniform overlap",
                 fontsize=13, fontweight="bold", y=1.00)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_dirichlet_99vs249.png",
                bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_dirichlet_99vs249.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print("== 3-way method comparison figures ==")
    fig_score_density()
    fig_rank_correlation()
    fig_score_volatility()
    fig_method_summary()
    fig_99vs249()
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
