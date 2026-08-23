#!/usr/bin/env python
"""Cross-method comparison figures for Dirichlet sensitivity analysis — enhanced.

Reads from projects/NeuralTF/results/:
  - dirichlet_top10_prioritized.csv          (centered)
  - dirichlet_uniform_top10.csv              (uniform)
  - dirichlet_uniform_full_rank.csv          (all 99 with both scores)
  - dirichlet_uniform_all249_full_rank.csv   (all 249 with uniform scores)
  - top10_neural_tfs_prioritized.csv         (fixed)

Outputs (PNGs into projects/NeuralTF/figures/):
  - fig_method_score_density.png          (KDE of all 99 candidates, 3 methods + stats)
  - fig_method_rank_correlation.png       (Spearman ρ heatmap + Pearson)
  - fig_method_score_volatility.png       (per-candidate score range + stats)
  - fig_method_summary.png                (5-panel: overlap, tracks, volatility, composite, concordance)
  - fig_method_99vs249.png                (99-neural vs 249-wide comparison)

Usage:
    python projects/NeuralTF/scripts/visualize_method_comparison.py
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
from scipy.stats import spearmanr, pearsonr, ks_2samp

REPO = Path(__file__).resolve().parents[3]
IN_DIR = REPO / "projects" / "NeuralTF" / "results"
OUT_DIR = REPO / "projects" / "NeuralTF" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Style
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

STREAMS = ["expression", "specificity", "reproducibility", "rnai",
           "correlation", "neural_enriched", "neural_specificity"]
W_DEFAULT = np.array([0.211, 0.105, 0.158, 0.158, 0.105, 0.158, 0.105])

C_FIXED = "#666666"
C_CENTERED = "#E69F00"
C_UNIFORM = "#009E73"
C_A = "#D55E00"
C_B = "#0072B2"
C_GRAY = "#999999"
C_GRID = "#E0E0E0"
C_BLACK = "#111111"
C_GREEN = "#009E73"
C_ORANGE = "#E69F00"
C_SKY = "#56B4E9"


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


def _build_all_scores() -> pd.DataFrame:
    """Return DataFrame with fixed, centered, uniform scores for all 99 candidates."""
    RUN = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
    rank = pd.read_csv(RUN / "rank_neural.csv")
    S = rank[STREAMS].to_numpy(dtype=float)
    mask = ~np.isnan(S)
    rank["fixed_score"] = np.where(mask, S, 0.0) @ W_DEFAULT

    full = pd.read_csv(IN_DIR / "dirichlet_uniform_full_rank.csv")
    full = full.rename(columns={"uniform_median_score": "uniform_score"})

    centered_csv = IN_DIR / "dirichlet_top10_prioritized.csv"
    if centered_csv.exists():
        centered = pd.read_csv(centered_csv)
        centered_map = dict(zip(centered["gene_id_v6"], centered["dirichlet_median_score"]))
        rank["centered_score"] = rank["gene_id"].map(centered_map)
    else:
        rank["centered_score"] = np.nan

    uniform_map = dict(zip(full["gene_id_v6"], full["uniform_score"]))
    rank["uniform_score"] = rank["gene_id"].map(uniform_map)

    out = rank[["gene_id", "gene_name", "fixed_score", "centered_score", "uniform_score"]].copy()
    return out.dropna()


def _add_stats_annotation(ax, data_groups, labels, x_pos=0.02, y_start=0.98):
    """Add statistical significance annotations between groups."""
    if len(data_groups) < 2:
        return
    y_max = max(max(g) for g in data_groups if len(g) > 0)
    y_step = 0.07 * (y_max if y_max > 0 else 1)
    y_pos = y_max + y_step
    from scipy.stats import mannwhitneyu, ks_2samp
    for i in range(len(data_groups)):
        for j in range(i + 1, len(data_groups)):
            if len(data_groups[i]) > 1 and len(data_groups[j]) > 1:
                try:
                    stat, p = mannwhitneyu(data_groups[i], data_groups[j], alternative='two-sided')
                    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                    ax.text(0.5, y_pos, f"{labels[i]} vs {labels[j]}: {sig} (p={p:.2e})",
                            ha='center', va='bottom', fontsize=7, color=C_BLACK, transform=ax.transAxes)
                    y_pos += 0.06
                except:
                    pass


# ---------------------------------------------------------------------------
# Figure 1: Score density KDE with statistics
# ---------------------------------------------------------------------------
def fig_score_density() -> None:
    df = _build_all_scores()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    _add_panel_letter(ax, "A")

    for col, color, label in [
        ("fixed_score", C_FIXED, "Fixed-weight"),
        ("centered_score", C_CENTERED, "Dirichlet-centered (k=40)"),
        ("uniform_score", C_UNIFORM, "Dirichlet-uniform (α=1)"),
    ]:
        vals = df[col].dropna()
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(vals, bw_method=0.2)
        x = np.linspace(vals.min() - 0.02, vals.max() + 0.02, 200)
        ax.plot(x, kde(x), color=color, linewidth=2, label=label)
        ax.fill_between(x, kde(x), color=color, alpha=0.15)
        # Add mean/median lines
        ax.axvline(vals.mean(), color=color, linestyle='--', linewidth=1, alpha=0.7)
        ax.axvline(vals.median(), color=color, linestyle=':', linewidth=1, alpha=0.7)

    ax.set_xlabel("Integrated score", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_title("A · Score density — 3 methods (99 neural candidates)",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, framealpha=0.9)
    _style_ax(ax)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_method_score_density.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_method_score_density.png")


# ---------------------------------------------------------------------------
# Figure 2: Spearman + Pearson correlation heatmap
# ---------------------------------------------------------------------------
def fig_rank_correlation() -> None:
    df = _build_all_scores()
    methods = ["fixed_score", "centered_score", "uniform_score"]
    labels = ["Fixed-weight", "Centered Dirichlet", "Uniform Dirichlet"]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    fig.subplots_adjust(wspace=0.3)

    # Panel A: Spearman
    ax = axes[0]
    _add_panel_letter(ax, "A")
    rho_spearman = np.zeros((3, 3))
    for i, ci in enumerate(["fixed_score", "centered_score", "uniform_score"]):
        for j, cj in enumerate(["fixed_score", "centered_score", "uniform_score"]):
            rho_spearman[i, j] = spearmanr(df[ci], df[cj]).correlation

    im1 = ax.imshow(rho_spearman, vmin=0.85, vmax=1.0, cmap="Reds", aspect="equal")
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{rho_spearman[i, j]:.3f}", ha="center", va="center",
                    fontsize=12, fontweight="bold",
                    color="white" if rho_spearman[i, j] > 0.95 else C_BLACK)
    ax.set_xticks(range(3))
    ax.set_yticks(range(3))
    ax.set_xticklabels(["Fixed", "Centered", "Uniform"], rotation=30, ha="right", fontsize=9)
    ax.set_yticklabels(["Fixed", "Centered", "Uniform"], fontsize=9)
    ax.set_title("A · Spearman ρ — rank correlation", fontsize=11, fontweight="bold")
    plt.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04, label="Spearman ρ")
    _style_ax(axes[0])

    # Panel B: Pearson
    ax = axes[1]
    _add_panel_letter(ax, "B")
    rho_pearson = np.zeros((3, 3))
    for i, ci in enumerate(["fixed_score", "centered_score", "uniform_score"]):
        for j, cj in enumerate(["fixed_score", "centered_score", "uniform_score"]):
            rho_pearson[i, j] = pearsonr(df[ci], df[cj])[0]

    im2 = ax.imshow(rho_pearson, vmin=0.85, vmax=1.0, cmap="Blues", aspect="equal")
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{rho_pearson[i, j]:.3f}", ha="center", va="center",
                    fontsize=12, fontweight="bold",
                    color="white" if rho_pearson[i, j] > 0.95 else C_BLACK)
    ax.set_xticks(range(3))
    ax.set_yticks(range(3))
    ax.set_xticklabels(["Fixed", "Centered", "Uniform"], rotation=30, ha="right", fontsize=9)
    ax.set_yticklabels(["Fixed", "Centered", "Uniform"], fontsize=9)
    ax.set_title("B · Pearson r — score correlation", fontsize=11, fontweight="bold")
    plt.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04, label="Pearson r")
    _style_ax(axes[1])

    fig.suptitle("Rank + score correlation across 3 methods (99 candidates)",
                 fontsize=12, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_method_rank_correlation.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_method_rank_correlation.png")


# ---------------------------------------------------------------------------
# Figure 3: Score volatility (per-candidate score range)
# ---------------------------------------------------------------------------
def fig_score_volatility() -> None:
    df = _build_all_scores()
    df["score_range"] = df[["fixed_score", "centered_score", "uniform_score"]].max(axis=1) - \
                        df[["fixed_score", "centered_score", "uniform_score"]].min(axis=1)
    df["score_cv"] = df[["fixed_score", "centered_score", "uniform_score"]].std(axis=1) / \
                     df[["fixed_score", "centered_score", "uniform_score"]].mean(axis=1)

    # Get top-20 by fixed score
    top20 = df.nlargest(20, "fixed_score").sort_values("score_range")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.subplots_adjust(wspace=0.3)

    # Panel A: Score trajectories
    ax = _add_panel_letter(_style_ax(plt.gca()), "A")
    plt.sca(_style_ax(plt.gca()))  # reset
    # Actually let's do this properly with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    _add_panel_letter(ax1, "A")

    centered_csv = IN_DIR / "dirichlet_top10_prioritized.csv"
    track_map = {}
    if (IN_DIR / "dirichlet_top10_prioritized.csv").exists():
        c = pd.read_csv(IN_DIR / "dirichlet_top10_prioritized.csv")
        track_map = dict(zip(c["gene_id_v6"], c["track"]))

    for i, (_, row) in enumerate(top20.iterrows()):
        gene_id = row["gene_id"]
        track = track_map.get(gene_id, "?")
        color = C_A if track == "A" else C_B
        scores = [row["fixed_score"], row["centered_score"], row["uniform_score"]]
        plt.sca(ax1)
        plt.plot(scores, [0, 1, 2], "o-", color=color, alpha=0.7,
                 linewidth=1.5, markersize=5, zorder=3)
        plt.text(-0.02, 0, f"{row['gene_name']}", ha="right", va="center",
                 fontsize=8, color=color, transform=plt.gca().get_yaxis_transform())
        plt.text(row["score_range"] + 0.005, 1, f"{row['score_range']:.3f}",
                 ha="left", va="center", fontsize=7, color=C_GRAY)

    ax1.set_xlim(df[["fixed_score", "centered_score", "uniform_score"]].min().min() - 0.02,
                 df[["fixed_score", "centered_score", "uniform_score"]].max().max() + 0.02)
    ax1.set_yticks([0, 1, 2])
    ax1.set_yticklabels(["Fixed", "Centered", "Uniform"])
    ax1.set_xlabel("Integrated score", fontsize=10)
    ax1.set_title("A · Score volatility — top 20 by fixed score",
                  fontsize=11, fontweight="bold")
    ax1.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_A,
               markersize=8, label="Track A (RNAi)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_B,
               markersize=8, label="Track B (novel)")
    ], loc="lower right", fontsize=8)
    _style_ax(ax1)

    # Panel B: Score range distribution
    _add_panel_letter(ax2, "B")
    all_ranges = df["score_range"].dropna()
    ax2.hist(all_ranges, bins=25, color=C_SKY, edgecolor=C_GRAY, alpha=0.75, density=True)
    ax2.axvline(all_ranges.mean(), color=C_ORANGE, linestyle='--', linewidth=1.5, 
                label=f'Mean = {all_ranges.mean():.4f}')
    ax2.axvline(all_ranges.median(), color=C_GREEN, linestyle=':', linewidth=1.5,
                label=f'Median = {all_ranges.median():.4f}')
    ax2.set_xlabel("Score range (max − min across 3 methods)", fontsize=10)
    ax2.set_ylabel("Density", fontsize=10)
    ax2.set_title("B · Score range distribution (all 99 candidates)",
                  fontsize=11, fontweight="bold")
    ax2.legend(fontsize=8, framealpha=0.9)
    _style_ax(ax2)

    fig.suptitle("Score volatility across 3 methods (99 neural candidates)",
                 fontsize=12, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_method_score_volatility.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_method_score_volatility.png")


# ---------------------------------------------------------------------------
# Figure 4: 5-panel summary (enhanced)
# ---------------------------------------------------------------------------
def fig_method_summary() -> None:
    """5-panel summary: overlap Venn, track assignment, volatility, composite, concordance."""
    fixed_csv = IN_DIR / "top10_neural_tfs_prioritized.csv"
    centered_csv = IN_DIR / "dirichlet_top10_prioritized.csv"
    uniform_csv = IN_DIR / "dirichlet_uniform_top10.csv"

    if not all(p.exists() for p in [fixed_csv, centered_csv, uniform_csv]):
        return

    fixed = pd.read_csv(fixed_csv)
    centered = pd.read_csv(centered_csv)
    uniform = pd.read_csv(uniform_csv)

    fixed_ids = set(fixed["gene_id_v6"])
    centered_ids = set(centered["gene_id_v6"])
    uniform_ids = set(uniform["gene_id_v6"])

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    fig.subplots_adjust(hspace=0.35, wspace=0.3)

    # Panel A: Venn diagram overlap
    ax = axes[0, 0]
    _add_panel_letter(ax, "A")
    from matplotlib_venn import venn3
    v = venn3([fixed_ids, centered_ids, uniform_ids],
              set_labels=("Fixed", "Centered", "Uniform"),
              set_colors=(C_FIXED, C_CENTERED, C_UNIFORM),
              alpha=0.6)
    for text in v.set_labels:
        if text: text.set_fontsize(10)
    for text in v.subset_labels:
        if text: text.set_fontsize(11)
    ax.set_title("A · Top-10 overlap across methods", fontsize=11, fontweight="bold")

    # Panel B: Track assignment consistency
    ax = axes[0, 1]
    _add_panel_letter(ax, "B")
    methods_data = [
        (fixed_ids, fixed, "Fixed"),
        (centered_ids, centered, "Centered"),
        (uniform_ids, uniform, "Uniform"),
    ]
    x = np.arange(3)
    width = 0.25
    for i, (ids, df, label) in enumerate(methods_data):
        a_count = len(df[df["track"] == "A"])
        b_count = len(df[df["track"] == "B"])
        ax.bar(x[i] - width, a_count, width, color=C_A, edgecolor="white",
               label="Track A" if i == 0 else "", zorder=3)
        ax.bar(x[i], b_count, width, color=C_B, edgecolor="white",
               label="Track B" if i == 0 else "", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(["Fixed", "Centered", "Uniform"], fontsize=9)
    ax.set_ylabel("Candidates in Top-10", fontsize=9)
    ax.set_title("B · Track assignment in Top-10", fontsize=11, fontweight="bold")
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc="upper right", fontsize=8)
    _style_ax(ax)

    # Panel C: Score volatility (top 20)
    ax = axes[1, 0]
    _add_panel_letter(ax, "C")
    df = _build_all_scores()
    df["score_range"] = df[["fixed_score", "centered_score", "uniform_score"]].max(axis=1) - \
                        df[["fixed_score", "centered_score", "uniform_score"]].min(axis=1)
    top20 = df.nlargest(20, "fixed_score").sort_values("score_range")

    centered_csv_path = IN_DIR / "dirichlet_top10_prioritized.csv"
    track_map = {}
    if centered_csv_path.exists():
        c = pd.read_csv(centered_csv_path)
        track_map = dict(zip(c["gene_id_v6"], c["track"]))

    for i, (_, row) in enumerate(top20.iterrows()):
        gene_id = row["gene_id"]
        track = track_map.get(gene_id, "?")
        color = C_A if track == "A" else C_B
        scores = [row["fixed_score"], row["centered_score"], row["uniform_score"]]
        ax.plot(scores, [0, 1, 2], "o-", color=color, alpha=0.7,
                linewidth=1.5, markersize=5, zorder=3)
        ax.text(-0.02, 0, f"{row['gene_name']}", ha="right", va="center",
                fontsize=8, color=color, transform=ax.get_yaxis_transform())
        ax.text(row["score_range"] + 0.005, 1, f"{row['score_range']:.3f}",
                ha="left", va="center", fontsize=7, color=C_GRAY)

    ax.set_xlim(df[["fixed_score", "centered_score", "uniform_score"]].min().min() - 0.02,
                df[["fixed_score", "centered_score", "uniform_score"]].max().max() + 0.02)
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["Fixed", "Centered", "Uniform"])
    ax.set_xlabel("Integrated score", fontsize=9)
    ax.set_title("C · Score volatility (top 20 by fixed score)",
                 fontsize=11, fontweight="bold")
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_A,
               markersize=8, label="Track A (RNAi)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_B,
               markersize=8, label="Track B (novel)")
    ], loc="lower right", fontsize=7)
    _style_ax(ax)

    # Panel D: Composite bonus effect (uniform top-10)
    ax = axes[1, 1]
    _add_panel_letter(ax, "D")
    uniform_top10 = pd.read_csv(IN_DIR / "dirichlet_uniform_top10.csv")
    fixed_top10 = pd.read_csv(IN_DIR / "top10_neural_tfs_prioritized.csv")

    for _, row in uniform_top10.iterrows():
        base = row["uniform_median_score"]
        comp = row["composite_score"]
        color = C_A if row["track"] == "A" else C_B
        marker = "o" if row["track"] == "A" else "^"
        ax.scatter(base, comp, color=color, s=60, marker=marker,
                   edgecolors="white", zorder=3)
        ax.annotate(row["gene_id_v6"].replace("dd_Smed_v6_", "dd"),
                    (base, comp), textcoords="offset points", xytext=(5, 3),
                    fontsize=7, color=C_BLACK)
    ax.plot([0, 1], [0, 1], "--", color=C_GRAY, linewidth=1.0, zorder=1,
            label="y = x (no bonus)")
    ax.set_xlim(0, 0.85)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("Base uniform-Dirichlet score", fontsize=9)
    ax.set_ylabel("Composite score (base + bonuses)", fontsize=9)
    ax.set_title("D · Composite bonus effect (uniform top-10)",
                 fontsize=11, fontweight="bold")
    ax.legend(handles=[
        Patch(facecolor=C_A, label="Track A (RNAi)"),
        Patch(facecolor=C_B, label="Track B (novel)"),
        Line2D([0], [0], color=C_GRAY, linewidth=1, linestyle="--", label="y = x (no bonus)"),
    ], loc="upper left", fontsize=8)
    _style_ax(ax)

    # Panel E: Three-way concordance (new panel)
    ax = axes[0, 2]
    _add_panel_letter(ax, "E")
    
    # Compute concordance: how often do methods agree on top-10 membership?
    all_top10_genes = fixed_ids | centered_ids | uniform_ids
    concordance_data = []
    for gene in all_top10_genes:
        in_fixed = gene in fixed_ids
        in_centered = gene in centered_ids
        in_uniform = gene in uniform_ids
        count = sum([in_fixed, in_centered, in_uniform])
        concordance_data.append({"gene": gene, "count": count})
    
    conc_df = pd.DataFrame(concordance_data)
    conc_counts = conc_df["count"].value_counts().sort_index()
    
    colors = [C_FIXED, C_CENTERED, C_UNIFORM]
    labels = ["Unique to 1 method", "Shared by 2 methods", "Consensus across all 3"]
    for i, (count, label, color) in enumerate(zip(conc_counts.values, labels, colors)):
        ax.bar(i, count, color=color, edgecolor="white", width=0.6, label=label, zorder=3)
        ax.text(i, count + 0.1, str(count), ha="center", va="bottom", fontsize=10, fontweight="bold")
    
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(labels, fontsize=9, rotation=15, ha="right")
    ax.set_ylabel("Number of genes", fontsize=9)
    ax.set_title("E · Top-10 consensus across methods", fontsize=11, fontweight="bold")
    ax.legend(fontsize=7, loc="upper right")
    _style_ax(ax)

    fig.suptitle("3-way method comparison summary — Dirichlet sensitivity analysis (99 candidates)",
                 fontsize=13, fontweight="bold", y=1.00)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_method_summary.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_method_summary.png")


# ---------------------------------------------------------------------------
# Figure 5: 99-neural vs 249-wide comparison (enhanced)
# ---------------------------------------------------------------------------
def fig_99vs249() -> None:
    """Compare 99-neural vs 249-wide uniform Dirichlet top-10."""
    unif99 = pd.read_csv(IN_DIR / "dirichlet_uniform_full_rank.csv")
    all249 = pd.read_csv(IN_DIR / "dirichlet_uniform_all249_full_rank.csv")
    unif99_top10 = pd.read_csv(IN_DIR / "dirichlet_uniform_top10.csv")
    unif99_top10_ids = set(unif99_top10["gene_id_v6"])
    all249_top10 = pd.read_csv(IN_DIR / "dirichlet_uniform_all249_top10.csv")
    all249_top10_ids = set(all249_top10["gene_id_v6"])

    RUN = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
    rank_neural = pd.read_csv(RUN / "rank_neural.csv")
    rank_all = pd.read_csv(RUN / "rank.csv")

    scores99 = dict(zip(rank_neural["gene_id"], unif99["uniform_median_score"]))
    scores249 = dict(zip(rank_all["gene_id"], all249["uniform_median_score"]))

    all99_ids = set(rank_neural["gene_id"])
    all249_ids_set = set(rank_all["gene_id"])

    rank99 = rank_neural.sort_values("integrated_score", ascending=False).reset_index(drop=True)
    rank99["rank"] = range(1, len(rank99) + 1)
    rank99_map = dict(zip(rank99["gene_id"], rank99["rank"]))

    rank249 = rank_all.sort_values("integrated_score", ascending=False).reset_index(drop=True)
    rank249["rank"] = range(1, len(rank249) + 1)
    rank249_map = dict(zip(rank249["gene_id"], rank249["rank"]))

    overlap_ids = all99_ids & all249_ids_set

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5.5), dpi=150)
    fig.subplots_adjust(wspace=0.35)

    # Panel A: Rank shift
    _add_panel_letter(ax_a, "A")
    rank_shift = []
    for gid in overlap_ids:
        r99 = rank99_map.get(gid, np.nan)
        r249 = rank249_map.get(gid, np.nan)
        if not np.isnan(r99) and not np.isnan(r249):
            rank_shift.append({"gene_id": gid, "rank_99": r99, "rank_249": r249,
                               "shift": r249 - r99})
    rank_shift_df = pd.DataFrame(rank_shift)

    colors = []
    for s in rank_shift_df["shift"]:
        if s > 0: colors.append(C_B)
        elif s < 0: colors.append(C_A)
        else: colors.append(C_GRAY)

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

    legend_elements = [
        Patch(facecolor=C_B, label="Rank higher (worse) in 249"),
        Patch(facecolor=C_A, label="Rank lower (better) in 249"),
        Patch(facecolor=C_GRAY, label="No change"),
    ]
    ax_a.legend(handles=legend_elements, loc="lower right", fontsize=7)

    # Panel B: Score comparison
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

    top10_99 = {g.replace("dd_Smed_v6_", "dd") for g in unif99_top10_ids}
    top10_249 = {g.replace("dd_Smed_v6_", "dd") for g in all249_top10_ids}

    for i, lbl in enumerate(labels):
        is_top10_99 = lbl in top10_99
        is_top10_249 = lbl in top10_249
        if is_top10_99 and is_top10_249:
            color = "#009E73"; marker = "o"; zord = 4
        elif is_top10_99:
            color = C_CENTERED; marker = "^"; zord = 3
        elif is_top10_249:
            color = C_A; marker = "s"; zord = 3
        else:
            color = C_GRAY; marker = "x"; zord = 2
        ax_b.scatter(x_scores[i], y_scores[i], color=color, s=50,
                     marker=marker, edgecolors="white", zorder=zord)
        if is_top10_99 or is_top10_249:
            ax_b.annotate(lbl, (x_scores[i], y_scores[i]),
                          textcoords="offset points", xytext=(5, 3),
                          fontsize=6, color=C_BLACK)

    # Add correlation
    rho, p_val = spearmanr(x_scores, y_scores)
    r_pearson, p_pearson = pearsonr(x_scores, y_scores)

    lims = [min(x_scores.min(), y_scores.min()) - 0.02,
            max(x_scores.max(), y_scores.max()) + 0.02]
    ax_b.plot(lims, lims, "--", color=C_GRAY, linewidth=1.0, zorder=1,
              label=f"y = x (ρ={rho:.3f}, r={r_pearson:.3f})")
    ax_b.set_xlim(lims)
    ax_b.set_ylim(lims)
    ax_b.set_xlabel("Uniform median score (99-neural scope)", fontsize=9)
    ax_b.set_ylabel("Uniform median score (249-wide scope)", fontsize=9)
    ax_b.set_title(f"B · Score comparison (ρ={rho:.3f}, r={r_pearson:.3f})",
                   fontsize=11, fontweight="bold")

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

    fig.suptitle("99-neural vs 249-wide Dirichlet-uniform overlap (uniform prior, α=1)",
                 fontsize=13, fontweight="bold", y=1.00)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_method_99vs249.png", bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_method_99vs249.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print("== 3-way method comparison figures (enhanced) ==")
    fig_score_density()
    fig_rank_correlation()
    fig_score_volatility()
    fig_method_summary()
    fig_99vs249()
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())