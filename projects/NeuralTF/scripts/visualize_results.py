#!/usr/bin/env python
"""Generate NeuralTF pipeline visualization figures from a run's outputs.

Reads only what the pipeline writes:
  - rank.csv (per-TF integrated scores + 7 evidence streams + proof_status + tier)
  - rank_neural.csv (neural-filtered subset)
  - figures/supplementary/go_gene_term_matrix_reduced.csv (GO dot plot)

Main figures (publication set):
  1. score_distributions  - per-stream score histograms + integrated
  2. candidate_summary    - 2x2: tiers, proof status, score by status, coverage
  3. top10_dual_track     - final Top-10, Track A (RNAi-validated) vs Track B (novel)
  4. evidence_heatmap     - top-30 evidence matrix (the core figure)
  5. candidate_funnel     - scored -> neural-filtered -> final candidates
  6. evidence_composition - stacked per-stream contribution of top-15
  7. stream_ablation      - rank sensitivity when each stream is removed
  8. top10_radar          - per-candidate 7-stream fingerprints
  9. go_dotplot           - GO-term coverage of the Top-10
10. integrated_vs_composite - composite (integrated + documented bonuses) vs
                                integrated for the final Top-10 (unique symbol
                                per candidate, track color, faint background)
  11. proof_status_violin - integrated-score distribution by proof status
  12. weight_sensitivity  - Top-10 rank bands under Dirichlet weight draws
  13. integrated_vs_neural_filter - ECDF of integrated score for all scored
                                  candidates, split by the neural filter

Figures are rendered in a Nature-publication style: minimal spines, explicit
informative titles, precise axis labels with units/definitions, Okabe-Ito
(color-vision-safe) palette, Arial sans-serif, 300 dpi PNGs.

Usage:
    python projects/NeuralTF/scripts/visualize_results.py
    python projects/NeuralTF/scripts/visualize_results.py \\
        --run projects/NeuralTF/runs/pipeline_run \\
        --out projects/NeuralTF/figures
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RUN = REPO_ROOT / "projects" / "NeuralTF" / "runs" / "pipeline_run"
DEFAULT_OUT = REPO_ROOT / "projects" / "NeuralTF" / "figures"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--run",
        type=Path,
        default=DEFAULT_RUN,
        help="Directory containing rank.csv / rank_neural.csv / pipeline_results.json",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUT,
        help="Directory to write figure PNGs to (created if missing)",
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# Nature-publication style
#
# - Arial (Nature body font) with DejaVu fallback for portability
# - minimal panels: left + bottom spines only, ticks outside
# - Okabe-Ito palette is colorblind-safe and used by Nature-style figures
# - 300 dpi exports, 3.5/7 in single/double-column width conventions
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
    "figure.dpi": 110,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "svg.fonttype": "none",  # keep text as text in vector exports
})

# Okabe-Ito (color-vision-safe)
C_BLUE = "#0072B2"
C_VERMILLION = "#D55E00"
C_GREEN = "#009E73"
C_SKY = "#56B4E9"
C_ORANGE = "#E69F00"
C_PURPLE = "#CC79A7"
C_GRAY = "#999999"
C_BLACK = "#000000"
C_YELLOW = "#F0E442"

TIER_COLORS = {"high": C_VERMILLION, "medium": C_BLUE, "low": C_GRAY}
PROOF_COLORS = {
    "known_rnai_validated": C_VERMILLION,
    "prior_fstf_not_tested": C_GREEN,
    "novel_candidate": C_BLUE,
}
STREAM_COLORS = {
    "expression": C_BLUE,
    "specificity": C_ORANGE,
    "reproducibility": C_GREEN,
    "rnai": C_VERMILLION,
    "correlation": C_PURPLE,
    "neural_enriched": C_SKY,
    "neural_specificity": C_GRAY,
}
SCORE_STREAMS = list(STREAM_COLORS.keys())

# Keep in sync with src/bioforge/evidence/scoring.py DEFAULT_WEIGHTS (the
# 0.05 weight freed by removing the FUNCTION stream was re-allotted
# proportionally, so all weights were scaled by 1/0.95 -> sum 1.0).
STREAM_WEIGHTS = {
    "expression": 0.211,
    "specificity": 0.105,
    "reproducibility": 0.158,
    "rnai": 0.158,
    "correlation": 0.105,
    "neural_enriched": 0.158,
    "neural_specificity": 0.105,
}

# Human-readable x-axis labels that name each scoring parameter precisely,
# so a reader understands the panel without consulting the methods section.
STREAM_XLABELS = {
    "expression": "Expression score\n(max log2 fold-change / 5, per atlas)",
    "specificity": "Specificity score\n(1 / number of supporting clusters)",
    "reproducibility": "Reproducibility score\n(atlases supporting / 3)",
    "rnai": "RNAi score\n(1 = target present in King 2024 mmc5 screen)",
    "correlation": "Correlation score\n(G0\u2013X1 TF-pair expression gain)",
    "neural_enriched": "Neural-enrichment score\n(1 = G0 subcluster log2FC \u2265 2)",
    "neural_specificity": "Neural-specificity score\n(1 / number of neural subclusters)",
    "integrated_score": "Integrated score (0\u20131,\nweighted sum of 7 streams)",
}

STREAM_SHORT = {
    "expression": "Expression",
    "specificity": "Specificity",
    "reproducibility": "Reproducibility",
    "rnai": "RNAi",
    "correlation": "Correlation",
    "neural_enriched": "Neural enriched",
    "neural_specificity": "Neural specificity",
}

# Stream order shared by the matrix figures (evidence heatmap + ablation):
# correlation is kept as the final column/row so the same stream ordering
# reads continuously across both heatmaps.
HEATMAP_STREAMS = [
    "expression",
    "specificity",
    "reproducibility",
    "rnai",
    "neural_enriched",
    "neural_specificity",
    "correlation",
]


def _tier_color(tier_name: str) -> str:
    return TIER_COLORS.get(str(tier_name).lower(), C_GRAY)


def _style_ax(ax) -> None:
    """Nature-style axis: only left/bottom spines, ticks facing out."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3.5, width=0.8)


def _panel_letter(ax, letter: str) -> None:
    """Bold italic panel letter in the top-left corner (Nature convention)."""
    ax.text(
        -0.18, 1.06, letter,
        transform=ax.transAxes, fontsize=12, fontweight="bold",
        fontstyle="italic", va="bottom", ha="left",
    )


def _save(fig, out_dir: Path, name: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"{name}.png"
    fig.tight_layout()
    fig.savefig(p, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {p.name}")


def _safe_get(df: pd.DataFrame, col: str) -> pd.Series:
    """Return df[col] if present, else a column of NaN with the right index."""
    if col in df.columns:
        return df[col]
    return pd.Series(np.nan, index=df.index, dtype=float)


def _label_for(df: pd.DataFrame) -> pd.Series:
    """Best identifier for plotting — gene_name if present else gene_id."""
    if "gene_name" in df.columns:
        return df["gene_name"].fillna(df.get("gene_id", pd.Series([""] * len(df)))).astype(str)
    return df.get("gene_id", pd.Series([""] * len(df))).astype(str)


def _renormalized(row: pd.Series, drop: str | None = None) -> float:
    """Integrated score over the streams present in `row`, optionally
    excluding `drop` — mirrors EvidenceScorer.integrated_score exactly."""
    num = den = 0.0
    for s in SCORE_STREAMS:
        if s == drop:
            continue
        v = row.get(s, np.nan)
        if pd.notna(v):
            num += STREAM_WEIGHTS[s] * v
            den += STREAM_WEIGHTS[s]
    return num / den if den > 0 else 0.0


def _final_shortlist(df_neural: pd.DataFrame | None) -> pd.DataFrame | None:
    """Merge the final dual-track shortlist (results/top10_neural_tfs_prioritized.csv,
    composite-based) with the neural rank table, so the publication Top-10
    figures show the *true final candidates* — not the integrated-score top-10,
    which can differ after annotation bonuses are applied."""
    csv = (REPO_ROOT / "projects" / "NeuralTF" / "results"
           / "top10_neural_tfs_prioritized.csv")
    if df_neural is None or df_neural.empty or not csv.exists():
        return None
    short = pd.read_csv(csv).rename(columns={"gene_id_v6": "gene_id"})
    if "gene_id" not in short.columns or "track" not in short.columns:
        return None
    out = short.merge(df_neural, on="gene_id", how="left", suffixes=("", "_r"))
    if out.empty:
        return None
    out["track"] = out["track"].astype(str).str.upper()
    order = {"A": 0, "B": 1}
    out["_track_ord"] = out["track"].map(order).fillna(2)
    out = out.sort_values(["_track_ord", "rank"]).drop(columns=["_track_ord"])
    return out


# ---------------------------------------------------------------------------
# Figure builders — return matplotlib Figure so the same function can be
# invoked by both this PNG-export script and the Streamlit UI (which can
# st.pyplot() a returned Figure directly without re-saving to disk).
# Signatures are load-bearing: the UI imports these by name.
# ---------------------------------------------------------------------------

def make_score_distributions(df: pd.DataFrame) -> plt.Figure | None:
    """Histograms of every evidence stream + the integrated score.

    Each panel names its scoring parameter on the x axis so the reader
    immediately knows what is plotted; the y axis is shared.
    """
    streams_present = [s for s in SCORE_STREAMS if s in df.columns]
    streams_present += ["integrated_score"]
    n = len(streams_present)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 3.7, nrows * 2.9))
    axes = np.atleast_1d(axes).ravel()
    letters = "abcdefghijklmnop"
    for i, ax in enumerate(axes):
        _style_ax(ax)
        if i >= n:
            ax.set_visible(False)
            continue
        s = streams_present[i]
        vals = df[s].dropna()
        if vals.empty:
            ax.set_visible(False)
            continue
        ax.hist(vals, bins=20, color=C_BLUE, alpha=0.85, edgecolor=C_BLACK,
                linewidth=0.5)
        ax.set_title(s.replace("_", " ").title(), fontsize=9.5)
        ax.set_xlabel(STREAM_XLABELS.get(s, "Score (0\u20131)"),
                      fontsize=7.5)
        ax.set_ylabel("Number of candidates")
        _panel_letter(ax, letters[i])
    fig.suptitle(
        f"Figure 1 \u2014 Distributions of the 7 evidence scores and the "
        f"integrated score across {len(df)} candidate TFs",
        fontsize=12, fontweight="bold", y=0.99,
    )
    fig.tight_layout()
    return fig


def make_candidate_summary(df: pd.DataFrame) -> plt.Figure | None:
    """2x2 summary: priority tiers, proof status, score by status, coverage."""
    has_tier = "tier" in df.columns
    has_proof = "proof_status" in df.columns
    if not (has_tier or has_proof):
        return None
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2))
    letters = ["a", "b", "c", "d"]

    # a) Priority tiers
    ax = axes[0, 0]
    _style_ax(ax)
    if has_tier:
        counts = df["tier"].value_counts()
        counts = counts.reindex(
            [t for t in ["high", "medium", "low"] if t in counts.index])
        colors = [_tier_color(t) for t in counts.index]
        ax.bar(counts.index, counts.values, color=colors, edgecolor="white",
               width=0.55)
        for i, v in enumerate(counts.values):
            ax.text(i, v + 1.5, str(v), ha="center", fontsize=8.5)
        ax.set_title("Priority tiers")
        ax.set_ylabel("Number of candidates")
    ax.set_xlabel("Tier")

    # b) Proof status
    ax = axes[0, 1]
    _style_ax(ax)
    if has_proof:
        counts = df["proof_status"].value_counts()
        counts = counts.reindex(
            [s for s in ["known_rnai_validated", "prior_fstf_not_tested",
                         "novel_candidate"] if s in counts.index])
        colors = [PROOF_COLORS.get(s, C_GRAY) for s in counts.index]
        ax.barh(counts.index, counts.values, color=colors, edgecolor="white",
                height=0.55)
        for i, v in enumerate(counts.values):
            ax.text(v + v * 0.02, i, str(v), va="center", fontsize=8.5)
        ax.set_xlim(0, max(counts.values) * 1.1)
        ax.set_title("Proof status")
        ax.set_xlabel("Number of candidates")
    ax.set_ylabel("Status")

    # c) Integrated score by proof status
    ax = axes[1, 0]
    _style_ax(ax)
    groups, labels = [], []
    if has_proof and "integrated_score" in df.columns:
        for s in ["known_rnai_validated", "prior_fstf_not_tested",
                  "novel_candidate"]:
            sub = df[df["proof_status"] == s]["integrated_score"].dropna()
            if not sub.empty:
                groups.append(sub.values)
                labels.append(s.replace("_", " ").title())
    if groups:
        ax.boxplot(groups, patch_artist=True, widths=0.5,
                   boxprops=dict(facecolor="#f0f0f0", edgecolor=C_BLACK,
                                 linewidth=0.8),
                   medianprops=dict(color=C_VERMILLION, linewidth=1.4),
                   whiskerprops=dict(color=C_BLACK, linewidth=0.8),
                   capprops=dict(color=C_BLACK, linewidth=0.8),
                   flierprops=dict(marker="o", markerfacecolor=C_BLUE,
                                   markeredgecolor="white", markersize=3))
        ax.set_xticklabels(labels, rotation=12, ha="right")
        ax.set_title("Integrated score by proof status")
        ax.set_ylabel("Integrated score (0\u20131)")
    else:
        ax.set_visible(False)

    # d) Evidence-stream coverage
    ax = axes[1, 1]
    _style_ax(ax)
    if "n_streams" in df.columns:
        counts = df["n_streams"].value_counts().reindex(
            range(1, 8), fill_value=0)
        ax.bar(counts.index, counts.values, color=C_BLUE, edgecolor="white",
               width=0.6)
        for x, v in counts.items():
            if v > 0:
                ax.text(x, v + 0.4, str(v), ha="center", fontsize=8)
        ax.set_xticks(range(1, 8))
        ax.set_title("Evidence-stream coverage")
        ax.set_xlabel("Supporting streams (of 7)")
        ax.set_ylabel("Number of candidates")
    else:
        ax.set_visible(False)

    for ax, letter in zip(axes.ravel(), letters):
        _panel_letter(ax, letter)
    fig.suptitle("Candidate summary \u2014 tiers, proof status, score "
                 "distribution, evidence coverage", fontsize=12,
                 fontweight="bold", y=0.99)
    fig.tight_layout()
    return fig


def make_top10_dual_track(df_neural: pd.DataFrame | None) -> plt.Figure | None:
    """Final Top-10: horizontal bars split into Track A (RNAi-validated)
    and Track B (novel). Uses the composite-based shortlist from
    prioritize_neural_tfs.py when available (final candidates), else the
    integrated-score top-10 of the neural subset as a fallback."""
    if df_neural is None or df_neural.empty:
        return None
    short = _final_shortlist(df_neural)
    if short is not None:
        top = short.iloc[::-1]
        score_col = "composite_score"
        score_label = ("Composite score (0\u20131, integrated + documented "
                       "annotation bonuses)")
        track_colors = {"A": PROOF_COLORS["known_rnai_validated"],
                        "B": PROOF_COLORS["novel_candidate"]}
        colors = [track_colors.get(str(t), C_GRAY) for t in top["track"]]
        title = ("Final Top-10 candidates \u2014 dual-track ranking\n"
                 "Track A (RNAi-validated) vs Track B (novel); N = "
                 f"{len(df_neural)} neural-filtered")
    else:
        top = (df_neural.sort_values("integrated_score", ascending=False)
               .head(10).iloc[::-1])
        score_col = "integrated_score"
        score_label = "Integrated score (0\u20131, weighted sum of 7 evidence streams)"
        colors = [PROOF_COLORS.get(str(s), C_GRAY)
                  for s in top.get("proof_status",
                                   pd.Series([""] * len(top)))]
        title = ("Top 10 neural candidates by integrated score (fallback "
                 "\u2014 run prioritize_neural_tfs.py for the final "
                 "composite shortlist)")
    labels = _label_for(top)
    fig, ax = plt.subplots(figsize=(7.5, max(4.6, 0.34 * len(top))))
    bars = ax.barh(labels, top[score_col], color=colors,
                   edgecolor="white", height=0.6)
    for bar, score in zip(bars, top[score_col]):
        ax.text(bar.get_width() + 0.008, bar.get_y() + bar.get_height() / 2,
                f"{score:.3f}", va="center", fontsize=7.5, color=C_BLACK)
    _style_ax(ax)
    ax.set_xlim(0, float(top[score_col].max()) * 1.18)
    ax.set_xlabel(score_label)
    ax.set_ylabel("Candidate gene")
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor=PROOF_COLORS["known_rnai_validated"], edgecolor="white",
              label="Track A \u2014 RNAi-validated"),
        Patch(facecolor=PROOF_COLORS["novel_candidate"], edgecolor="white",
              label="Track B \u2014 novel candidate"),
    ]
    if short is None and (top.get("proof_status") == "prior_fstf_not_tested").any():
        handles.append(Patch(facecolor=PROOF_COLORS["prior_fstf_not_tested"],
                             edgecolor="white",
                             label="Prior FSTF (literature)"))
    ax.legend(handles=handles, loc="lower right", frameon=True,
              framealpha=0.9, edgecolor="0.8")
    ax.set_title(title, fontsize=10.5)
    fig.tight_layout()
    return fig


def make_evidence_heatmap(df: pd.DataFrame, n: int = 30) -> plt.Figure | None:
    """Heatmap of evidence-stream scores for the top-N candidates."""
    streams = [s for s in HEATMAP_STREAMS if s in df.columns]
    if not streams:
        return None
    top = df.head(n)
    heat = top.set_index(_label_for(top))[streams]
    fig, ax = plt.subplots(figsize=(7.5, max(4.8, 0.3 * len(heat))))
    sns.heatmap(
        heat.astype(float), annot=True, fmt=".2f", cmap="RdYlGn", center=0.5,
        annot_kws={"fontsize": 6.5},
        cbar_kws={"label": "Stream score (0\u20131)",
                  "shrink": 0.6, "aspect": 24},
        linewidths=0.4, linecolor="white", ax=ax,
    )
    ax.set_title(
        f"Evidence-stream scores of the top {min(n, len(heat))} candidate TFs\n"
        "(rows = candidates, columns = 7 weighted evidence streams)",
        fontsize=10.5,
    )
    ax.set_xlabel("Evidence stream")
    ax.set_ylabel("Candidate gene")
    ax.tick_params(axis="both", length=0)
    fig.tight_layout()
    return fig


def make_evidence_heatmap_all_genes(df: pd.DataFrame) -> plt.Figure | None:
    """Heatmap of evidence-stream scores for ALL neural candidates (99 genes)."""
    streams = [s for s in HEATMAP_STREAMS if s in df.columns]
    if not streams:
        return None
    # Sort by integrated_score descending
    df_sorted = df.sort_values("integrated_score", ascending=False).reset_index(drop=True)
    heat = df_sorted.set_index(_label_for(df_sorted))[streams]

    n_genes = len(heat)
    fig_height = max(14, 0.28 * n_genes)
    fig, ax = plt.subplots(figsize=(9, fig_height))

    # Custom colormap: white -> light green -> dark green
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "nature_green", ["#FFFFFF", "#E8F5E9", "#2E7D32"])

    im = ax.imshow(heat.astype(float).to_numpy(), aspect="auto", cmap=cmap,
                   vmin=0, vmax=1, interpolation="nearest")

    # Y-axis: gene labels
    ax.set_yticks(np.arange(n_genes))
    ax.set_yticklabels(heat.index.tolist(), fontsize=5.5, fontfamily="sans-serif")
    ax.set_ylabel("Candidate gene (ranked by integrated score)", fontsize=10, labelpad=8)

    # X-axis: stream labels using STREAM_SHORT
    ax.set_xticks(np.arange(len(streams)))
    ax.set_xticklabels([STREAM_SHORT[s] for s in streams],
                       fontsize=8, rotation=0, ha="center", fontfamily="sans-serif")
    ax.set_xlabel("Evidence stream", fontsize=10, labelpad=8)

    # Colorbar on the right, outside plot
    cbar = fig.colorbar(im, ax=ax, shrink=0.5, aspect=30, pad=0.02,
                        location="right")
    cbar.set_label("Stream score (0\u20131)", fontsize=9, labelpad=6)
    cbar.ax.tick_params(labelsize=8)

    # White grid lines
    ax.set_xticks(np.arange(len(streams) + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(n_genes + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.5)
    ax.tick_params(which="minor", length=0)

    # Annotate top 20 with values
    top_annotate = min(20, n_genes)
    heat_vals = heat.astype(float).to_numpy()
    for i in range(top_annotate):
        for j in range(len(streams)):
            val = heat_vals[i, j]
            if not np.isnan(val):
                text_color = "white" if val > 0.6 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=4.5, color=text_color, fontweight="medium")

    # Title
    ax.set_title("Evidence-stream scores for all 99 neural-filtered candidates\n"
                 "(rows = candidates sorted by integrated score; columns = 7 evidence streams)",
                 fontsize=11, fontweight="bold", pad=14)

    fig.tight_layout(rect=[0, 0, 0.85, 1])  # Leave space for colorbar
    return fig


def make_candidate_funnel(df: pd.DataFrame,
                          df_neural: pd.DataFrame | None) -> plt.Figure | None:
    """Funnel: scored candidates -> neural-filtered -> final Top-10."""
    n_all = len(df)
    n_neural = len(df_neural) if df_neural is not None else 0
    short = _final_shortlist(df_neural)
    n_top = len(short) if short is not None else min(10, n_neural)
    stages = [
        "Detectably enriched TFs scored",
        "Neural-filtered candidates\n(neural-G0 log2FC \u2265 2 in \u22651 of 77 "
        "neural subclusters)",
        "Final candidates\n(top 5 per track)",
    ]
    values = [n_all, n_neural, n_top]
    widths = np.array(values, dtype=float) / max(values) * 8.0
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for i in range(3):
        left = (8.0 - widths[i]) / 2.0
        ax.barh(2 - i, widths[i], left=left, height=0.62,
                color=[C_GRAY, C_BLUE, C_VERMILLION][i],
                edgecolor="white")
        ax.text(4.0, 2 - i, f"{values[i]:,}", ha="center", va="center",
                fontsize=11, fontweight="bold", color="white")
    for i in range(2):
        pct = 100.0 * values[i + 1] / values[i]
        ax.text(4.0, 1.28 - i * 1.0, f"\u2193 {pct:.1f}% retained",
                ha="center", fontsize=8, color="#444444")
    _style_ax(ax)
    ax.set_yticks([2, 1, 0])
    ax.set_yticklabels(stages, fontsize=9)
    ax.set_xlim(0, 8)
    ax.set_xticks([])
    ax.set_title("Candidate funnel \u2014 from scored TFs to final candidates",
                 fontsize=11)
    fig.tight_layout()
    return fig


def make_evidence_composition(df: pd.DataFrame, n: int = 15) -> plt.Figure | None:
    """Stacked horizontal bars: which streams drive each candidate's score.

    Segments are the *weighted* contributions w(s) x s, so the stacked
    total is the un-normalized weighted sum (0-1 scale) — for candidates
    with all 7 streams present it equals the integrated score exactly.
    """
    streams = [s for s in SCORE_STREAMS
               if s in df.columns and df[s].notna().any()]
    if not streams:
        return None
    top = df.head(n).iloc[::-1]
    labels = _label_for(top)
    fig, ax = plt.subplots(figsize=(7.5, max(4.2, 0.38 * len(top))))
    bottom = np.zeros(len(top))
    for s in streams:
        vals = top[s].fillna(0.0).to_numpy(dtype=float) \
            * STREAM_WEIGHTS.get(s, 1.0)
        ax.barh(labels, vals, left=bottom, label=s,
                color=STREAM_COLORS[s], edgecolor="white", linewidth=0.3,
                height=0.6)
        bottom += vals
    _style_ax(ax)
    ax.set_xlim(0, max(bottom.max() * 1.15, 1.0))
    ax.set_xlabel("Weighted evidence contribution (0\u20131,\n"
                  "w\u00b7score per stream; total = un-normalized weighted sum)")
    ax.set_ylabel("Candidate gene")
    ax.set_title(
        f"Evidence-stream composition of the top {min(n, len(top))} candidates",
        fontsize=10.5,
    )
    ax.legend(loc="lower right", frameon=True, framealpha=0.9,
              edgecolor="0.8", ncol=2, fontsize=7.5)
    fig.tight_layout()
    return fig


def make_stream_ablation(df: pd.DataFrame, n: int = 30) -> plt.Figure | None:
    """Rank change per candidate when each stream is removed from the score.

    For every stream the integrated score is recomputed over the remaining
    streams (renormalized weights, identical to EvidenceScorer), the ranking
    is re-derived, and Delta = rank_without_stream - rank_full is plotted.
    Red = candidate ranks lower without the stream (the stream was
    load-bearing for it); blue = ranks higher; white = unchanged.
    """
    if "integrated_score" not in df.columns:
        return None
    rank_full = df["integrated_score"].rank(ascending=False, method="min")
    top_idx = df.nlargest(n, "integrated_score").index
    rows = {}
    for s in HEATMAP_STREAMS:
        ablated = df.apply(lambda r: _renormalized(r, drop=s), axis=1)
        rank_ab = ablated.rank(ascending=False, method="min")
        rows[s] = (rank_ab - rank_full).loc[top_idx].values
    heat = pd.DataFrame(rows, index=_label_for(df.loc[top_idx])).T
    heat.index = [STREAM_SHORT[s] for s in HEATMAP_STREAMS]
    vmax = max(float(heat.abs().values.max()), 1.0)
    fig, ax = plt.subplots(figsize=(12.5, 4.0))
    sns.heatmap(
        heat.astype(float), annot=True, fmt=".0f", cmap="RdBu_r", center=0,
        vmin=-vmax, vmax=vmax, linewidths=0.4, linecolor="white",
        annot_kws={"fontsize": 6},
        cbar_kws={"label": "\u0394rank (rank without stream \u2212 full rank)",
                  "shrink": 0.7, "aspect": 22},
        ax=ax,
    )
    ax.set_title(
        f"Rank sensitivity to each evidence stream \u2014 top {min(n, len(top_idx))} "
        "candidates\nred = candidate ranks lower without the stream; blue = higher; "
        "white = no change",
        fontsize=10.5,
    )
    ax.set_xlabel("Candidate gene (ordered by full integrated score)")
    ax.set_ylabel("Removed evidence stream")
    ax.tick_params(axis="both", length=0)
    fig.tight_layout()
    return fig


def make_top10_radar(df_neural: pd.DataFrame | None) -> plt.Figure | None:
    """Radar fingerprints of the final Top-10 candidates across the 7 streams."""
    if df_neural is None or df_neural.empty:
        return None
    short = _final_shortlist(df_neural)
    top = short if short is not None else (
        df_neural.sort_values("integrated_score", ascending=False).head(10))
    if top.empty:
        return None
    n_streams = len(SCORE_STREAMS)
    angles = np.linspace(0, 2 * np.pi, n_streams, endpoint=False).tolist()
    angles_closed = angles + angles[:1]
    labels = [STREAM_SHORT[s] for s in SCORE_STREAMS]
    fig, axes = plt.subplots(
        2, 5, subplot_kw=dict(polar=True), figsize=(13.5, 6.6))
    for i, (_, row) in enumerate(top.iterrows()):
        ax = axes.ravel()[i]
        vals = [0.0 if pd.isna(row.get(s)) else float(row[s])
                for s in SCORE_STREAMS]
        vals_closed = vals + vals[:1]
        ax.plot(angles_closed, vals_closed, color=C_BLUE, linewidth=1.1)
        ax.fill(angles_closed, vals_closed, color=C_BLUE, alpha=0.16)
        ax.set_xticks(angles)
        ax.set_xticklabels(labels, fontsize=5.6)
        ax.set_ylim(0, 1)
        ax.set_yticks([0.5])
        ax.set_yticklabels([])
        ax.grid(True, color="0.85", linewidth=0.5)
        score = row.get("composite_score", row.get("integrated_score", float("nan")))
        ax.set_title(
            f"{_label_for(pd.DataFrame([row])).iloc[0]}  "
            f"{score:.3f}",
            fontsize=8, pad=10,
        )
    fig.suptitle(
        "Final Top-10 candidate evidence fingerprints (7 streams; title = gene, "
        "composite score)",
        fontsize=12, fontweight="bold", y=1.0,
    )
    fig.tight_layout()
    return fig


def make_go_dotplot(df_neural: pd.DataFrame | None) -> plt.Figure | None:
    """Dot plot of GO-term coverage among the final Top-10 candidates.

    Reads figures/supplementary/go_gene_term_matrix_reduced.csv (produced by
    make_supp_go_figures.py); returns None when the matrix is unavailable.
    """
    matrix_csv = (REPO_ROOT / "projects" / "NeuralTF" / "figures"
                  / "supplementary" / "go_gene_term_matrix_reduced.csv")
    if df_neural is None or df_neural.empty or not matrix_csv.exists():
        return None
    short = _final_shortlist(df_neural)
    top = short if short is not None else (
        df_neural.sort_values("integrated_score", ascending=False).head(10))
    names = _label_for(top).tolist()
    top_ids = top["gene_id"].astype(str).tolist()
    go = pd.read_csv(matrix_csv)
    keep = go[go["gene_id"].astype(str).isin(set(top_ids))]
    if keep.empty:
        return None
    terms = [c for c in keep.columns if c.startswith("GO:")]
    if not terms:
        return None
    freq = keep[terms].sum()
    freq = freq[freq >= 2].sort_values(ascending=False).head(12)
    if freq.empty:
        freq = keep[terms].sum().sort_values(ascending=False).head(12)
    rows = list(freq.index)
    fig, ax = plt.subplots(figsize=(9.5, max(4.5, 0.42 * len(rows))))
    ref_csv = REPO_ROOT / "projects" / "NeuralTF" / "figures" / "go_term_reference.csv"
    cat_of: dict[str, str] = {}
    if ref_csv.exists():
        ref = pd.read_csv(ref_csv)
        for _, r in ref.iterrows():
            n_, t_ = str(r.get("neural_go", "")), str(r.get("tf_go", ""))
            if n_ == "yes" and t_ == "yes":
                c = "both"
            elif n_ == "yes":
                c = "neural"
            elif t_ == "yes":
                c = "tf"
            else:
                c = "other"
            cat_of[str(r.get("go_id", ""))] = c
    SHAPE = {"neural": "o", "tf": "s", "both": "D", "other": "v"}
    SHAPE_LABEL = {"neural": "neural GO term",
                   "tf": "TF-related GO term",
                   "both": "neural + TF GO term",
                   "other": "other GO term"}
    cmap = plt.cm.viridis
    norm = plt.Normalize(vmin=float(freq.min()), vmax=float(freq.max()))
    for j, term in enumerate(rows):
        marker = SHAPE.get(cat_of.get(term, "other"), "o")
        for i, gid in enumerate(top_ids):
            row = keep[keep["gene_id"].astype(str) == gid]
            if row.empty or row.iloc[0].get(term, 0) != 1:
                continue
            ax.scatter(i, j, s=30 + 45 * int(freq[term]),
                       c=[cmap(norm(int(freq[term])))],
                       marker=marker, alpha=0.9,
                       edgecolor="white", linewidth=0.4)
    _style_ax(ax)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(
        [f"{t} ({int(freq[t])})" for t in rows], fontsize=7.5)
    ax.set_xlabel("Candidate gene (final Top-10 ordering)")
    ax.set_ylabel("GO term (count of Top-10 sharing it)")
    ax.set_xlim(-0.6, len(names) - 0.4)
    from matplotlib.lines import Line2D
    handles = []
    for cat, m in SHAPE.items():
        handles.append(Line2D([0], [0], marker=m, linestyle="",
                              markersize=8, markerfacecolor="0.75",
                              markeredgecolor="black",
                              label=SHAPE_LABEL[cat]))
    ax.legend(handles=handles, loc="lower left", frameon=True,
              framealpha=0.9, edgecolor="0.8", fontsize=7.5)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("Top-10 genes sharing the term", fontsize=8)
    ax.set_title(
        "GO-term coverage of the Top-10 neural candidates\n"
        "(marker shape = GO-term category; marker color = how many of the "
        "Top-10 share it)",
        fontsize=10.5,
    )
    fig.tight_layout()
    return fig


def make_integrated_vs_composite(df_neural: pd.DataFrame | None) -> plt.Figure | None:
    """Top-10 scatter of composite vs integrated score.

    composite = integrated + documented annotation bonuses (TF domain +0.05,
    neural GO +0.03, TF GO +0.02, human ortholog +0.02, RNAi +0.02). Reuses
    the exact bonus logic from bioforge/projects/neuraltf/prioritize.py on
    the PlanMine parquet + King mmc4 catalog. The final Top-10 is drawn with
    a unique symbol per candidate (color = track); the other neural
    candidates are faint background.
    """
    if df_neural is None or df_neural.empty:
        return None
    try:
        from bioforge.projects.neuraltf.prioritize import (
            compute_composite, merge_annotations, prepare_candidates,
            summarize_annotations,
        )
        king_dir = REPO_ROOT / "datasets" / "raw" / "Supplementary_Data_ King_2024"
        mmc4_csv = sorted(king_dir.glob("*mmc4*.xlsx"))[0]
        mmc4 = pd.read_excel(mmc4_csv, sheet_name="TF")
        parquet = REPO_ROOT / "datasets" / "processed" / "planmine_annotations.parquet"
        if not parquet.exists():
            return None
        ann = summarize_annotations(pd.read_parquet(parquet))
        comp = compute_composite(
            merge_annotations(prepare_candidates(df_neural, mmc4), ann))
    except Exception:
        return None
    short = _final_shortlist(df_neural)
    short_ids = set(short["gene_id"].astype(str)) if short is not None else set()
    comp = comp.reset_index(drop=True)
    x = pd.to_numeric(comp["integrated_score"], errors="coerce")
    y = pd.to_numeric(comp["composite_score"], errors="coerce")
    ok = x.notna() & y.notna()
    top10_csv = (REPO_ROOT / "projects" / "NeuralTF" / "results"
                 / "top10_neural_tfs_prioritized.csv")
    order: list[str] = []
    track_of: dict[str, str] = {}
    if top10_csv.exists():
        t10 = pd.read_csv(top10_csv).sort_values("rank")
        order = [str(g) for g in t10["gene_id_v6"]]
        track_of = dict(zip(order, [str(t) for t in t10["track"]]))
    if not order:
        order = sorted(short_ids)
    MARKERS = ["o", "v", "^", "<", ">", "s", "p", "D", "P", "X"]
    TRACK_COLOR = {"A": C_ORANGE, "B": C_BLUE}

    fig, ax = plt.subplots(figsize=(7.6, 5.9))
    bg = comp[ok & ~comp["gene_id"].astype(str).isin(short_ids)]
    if not bg.empty:
        ax.scatter(pd.to_numeric(bg["integrated_score"], errors="coerce"),
                   pd.to_numeric(bg["composite_score"], errors="coerce"),
                   s=14, color=C_GRAY, alpha=0.3, edgecolor="none",
                   zorder=1, label=f"other neural candidates (n={len(bg)})")
    top_by_id = {str(r["gene_id"]): r for _, r in comp.iterrows()
                 if str(r["gene_id"]) in short_ids}
    for i, gid in enumerate(order):
        if gid not in top_by_id:
            continue
        row = top_by_id[gid]
        marker = MARKERS[i % len(MARKERS)]
        color = TRACK_COLOR.get(track_of.get(gid, "B"), C_BLUE)
        ax.scatter(pd.to_numeric(row["integrated_score"], errors="coerce"),
                   pd.to_numeric(row["composite_score"], errors="coerce"),
                   s=85, marker=marker, color=color,
                   edgecolor=C_BLACK, linewidth=0.8, zorder=5)
        ax.annotate(row["gene_name"],
                    (row["integrated_score"], row["composite_score"]),
                    textcoords="offset points", xytext=(7, 5),
                    fontsize=7, color=C_BLACK, zorder=6)
    lims = (0.0, max(x.max(), y.max(), 1.0) * 1.02)
    ax.plot(lims, lims, color=C_GRAY, linestyle="--", linewidth=0.9,
            zorder=1, label="y = x (no bonuses)")
    _style_ax(ax)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Integrated score (0\u20131, 7 weighted evidence streams)")
    ax.set_ylabel("Composite score (0\u20131, integrated + documented bonuses)")
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=6,
                      markerfacecolor=C_GRAY, markeredgecolor="none",
                      label=f"other neural candidates (n={len(bg)})")]
    for i, gid in enumerate(order):
        if gid not in top_by_id:
            continue
        row = top_by_id[gid]
        trk = track_of.get(gid, "B")
        handles.append(Line2D(
            [0], [0], marker=MARKERS[i % len(MARKERS)], linestyle="",
            markersize=8, markerfacecolor=TRACK_COLOR.get(trk, C_BLUE),
            markeredgecolor=C_BLACK, label=f"{row['gene_name']} (Track {trk})"))
    ax.legend(handles=handles, loc="upper left", frameon=True,
              framealpha=0.92, edgecolor="0.8", fontsize=7)
    ax.set_title(
        "Annotation bonuses for the final Top-10\n"
        "symbol = candidate (color = track); faint = other neural candidates",
        fontsize=10.5,
    )
    fig.tight_layout()
    return fig


def make_neural_filter_vs_score(df: pd.DataFrame) -> plt.Figure | None:
    """Overlapping proportion histograms of integrated score for every scored
    candidate, split by the neural filter (neural-enriched G0 signal OR RNAi
    validation).

    Answers: does the neural filter select candidates with higher
    integrated scores?
    """
    need = {"integrated_score", "neural_enriched", "rnai"}
    if df is None or df.empty or not need.issubset(df.columns):
        return None
    score = pd.to_numeric(df["integrated_score"], errors="coerce")
    rnai = pd.to_numeric(df["rnai"], errors="coerce").fillna(0)
    neural = df["neural_enriched"].notna() | (rnai > 0)
    if int(neural.sum()) == 0 or int((~neural).sum()) == 0:
        return None
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    bins = np.linspace(0.0, 1.0, 26)
    stats: dict[str, float] = {}
    for mask, color, label in [
            (~neural, C_GRAY, "Excluded (not neural-filtered)"),
            (neural, C_BLUE, "Neural-filtered")]:
        vals = score[mask].dropna().to_numpy()
        if vals.size == 0:
            continue
        ax.hist(vals, bins=bins, density=False, alpha=0.55, color=color,
                edgecolor="white", linewidth=0.5, zorder=2,
                label=f"{label} (n={vals.size})")
        med = float(np.median(vals))
        stats[label] = med
        ax.axvline(med, color=color, linestyle="--", lw=1.1, alpha=0.9,
                   zorder=3)
    _style_ax(ax)
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Integrated score (0\u20131, 7 weighted evidence streams)")
    ax.set_ylabel("Number of candidates")
    if len(stats) == 2:
        ax.text(0.99, 0.97,
                f"median \u2014 neural: {stats.get('Neural-filtered', float('nan')):.3f}  "
                f"excluded: {stats.get('Excluded (not neural-filtered)', float('nan')):.3f}",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=8, color="0.25")
    ax.legend(loc="upper right", frameon=True, framealpha=0.92,
              edgecolor="0.8", fontsize=8)
    ax.set_title(
        "Integrated score vs neural filter (all scored candidates)\n"
        "neural filter = neural-enriched G0 signal or RNAi validation",
        fontsize=10.5,
    )
    fig.tight_layout()
    return fig


def make_proof_status_violin(df: pd.DataFrame) -> plt.Figure | None:
    """Integrated-score distribution by proof status (violin + jitter)."""
    if df is None or df.empty or "integrated_score" not in df.columns \
            or "proof_status" not in df.columns:
        return None
    order = [s for s in ["known_rnai_validated", "prior_fstf_not_tested",
                         "novel_candidate"]
             if s in df["proof_status"].unique()]
    if not order:
        return None
    labels = [s.replace("_", " ").title() + f"\n(n = {len(df[df['proof_status'] == s])})"
              for s in order]
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    parts = ax.violinplot(
        [df[df["proof_status"] == s]["integrated_score"].dropna().values
         for s in order],
        positions=range(1, len(order) + 1), widths=0.7, showmeans=False,
        showmedians=False, showextrema=False,
    )
    for body, s in zip(parts["bodies"], order):
        body.set_facecolor(PROOF_COLORS.get(s, C_GRAY))
        body.set_alpha(0.45)
        body.set_edgecolor(PROOF_COLORS.get(s, C_GRAY))
    rng = np.random.default_rng(7)
    for i, s in enumerate(order, start=1):
        vals = df[df["proof_status"] == s]["integrated_score"].dropna().values
        jitter = rng.normal(0, 0.06, size=len(vals))
        ax.scatter(i + jitter, vals, s=14, color=PROOF_COLORS.get(s, C_GRAY),
                   alpha=0.65, edgecolor="white", linewidth=0.3, zorder=3)
        ax.axhline(np.median(vals), color=PROOF_COLORS.get(s, C_GRAY),
                   linestyle="--", linewidth=1.0, alpha=0.7)
    _style_ax(ax)
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Integrated score (0\u20131)")
    ax.set_title(
        "Integrated-score distribution by proof status\n"
        "validated TFs carry the RNAi signal; novel candidates are the "
        "untested tail",
        fontsize=10.5,
    )
    fig.tight_layout()
    return fig


def make_weight_sensitivity(df_neural: pd.DataFrame | None,
                            out_dir: Path | None = None) -> plt.Figure | None:
    """Top-10 rank bands when stream weights are perturbed.

    For every draw a fresh weight vector is sampled from a Dirichlet
    distribution concentrated on the 7 default weights (alpha = weight x 40,
    absent streams zeroed — mirrors the scorer's renormalization), all
    neural candidates are re-scored and re-ranked, and the ranks of the
    final Top-10 members are recorded. Whiskers = min-max rank over 1000
    draws; the dashed line marks the top-10 cutoff.

    When ``out_dir`` is given, the full draw data is persisted next to the
    figure:

    - ``weight_sensitivity_draws.csv`` - long table (draw, gene_id,
      gene_name, rank, in_top_10) for every candidate x draw;
    - ``weight_sensitivity_top10_challengers.csv`` - per-candidate summary
      marking who replaced the baseline Top-10 (counts + which baseline
      member was displaced most often).
    """
    if df_neural is None or df_neural.empty:
        return None
    short = _final_shortlist(df_neural)
    if short is None or len(short) < 2:
        return None
    streams = [s for s in SCORE_STREAMS if s in df_neural.columns]
    S = df_neural[streams].to_numpy(dtype=float)
    mask = ~np.isnan(S)
    W = np.array([STREAM_WEIGHTS[s] for s in streams], dtype=float)
    rng = np.random.default_rng(2024)
    n_draws = 1000
    k = 40.0
    n = len(df_neural)
    ranks = np.empty((n_draws, n), dtype=int)
    for d in range(n_draws):
        alpha = W[None, :] * mask * k + 1e-9
        w = rng.gamma(alpha, 1.0)
        w = w / w.sum(axis=1, keepdims=True)
        score = np.nansum(w * np.where(mask, S, 0.0), axis=1)
        order = np.argsort(-score, kind="stable")
        ranks[d, order] = np.arange(1, n + 1)
    if out_dir is not None:
        _persist_weight_draws(df_neural, ranks, short, out_dir)
    idx_of = {gid: i for i, gid in
              enumerate(df_neural["gene_id"].astype(str))}
    short = short.reset_index(drop=True)
    rows = []
    for _, row in short.iterrows():
        i = idx_of.get(str(row.get("gene_id")), None)
        if i is None:
            continue
        r = ranks[:, i]
        rows.append({
            "label": _label_for(pd.DataFrame([row])).iloc[0],
            "track": str(row.get("track", "")).upper(),
            "median": float(np.median(r)),
            "lo": float(r.min()),
            "hi": float(r.max()),
            "drops_out": bool(r.max() > 10),
        })
    if not rows:
        return None
    d = pd.DataFrame(rows)
    track_colors = {"A": PROOF_COLORS["known_rnai_validated"],
                    "B": PROOF_COLORS["novel_candidate"]}
    colors = [track_colors.get(t, C_GRAY) for t in d["track"]]
    x = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    rng = np.random.default_rng(11)
    for i, (xi, r) in enumerate(zip(x, rows)):
        ax.vlines(xi, r["lo"], r["hi"], color=colors[i], linewidth=2.2,
                  alpha=0.85)
        draws = ranks[:, idx_of[str(short.loc[i, "gene_id"])]]
        ax.scatter(xi + rng.normal(0, 0.06, size=len(draws)), draws,
                   s=4, alpha=0.22, color=colors[i], zorder=2)
        ax.scatter(xi, r["median"], s=30, color="white",
                   edgecolor=colors[i], linewidth=1.4, zorder=4)
    ax.axhline(10, color=C_GRAY, linestyle="--", linewidth=1.0,
               label="top-10 cutoff")
    _style_ax(ax)
    ax.set_xticks(x)
    ax.set_xticklabels(d["label"], rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("Rank among all neural candidates (lower = better)")
    ax.set_ylim(len(df_neural) + 1, 0)
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor=PROOF_COLORS["known_rnai_validated"], edgecolor="white",
              label="Track A \u2014 RNAi-validated"),
        Patch(facecolor=PROOF_COLORS["novel_candidate"], edgecolor="white",
              label="Track B \u2014 novel candidate"),
        plt.Line2D([0], [0], color=C_GRAY, linestyle="--",
                   label="top-10 cutoff"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=True,
              framealpha=0.9, edgecolor="0.8", fontsize=7.5)
    ax.set_title(
        "Top-10 rank stability under weight uncertainty\n"
        "(1000 Dirichlet draws around the 7 default weights; whisker = "
        "min\u2013max rank, dot = median)",
        fontsize=10.5,
    )
    fig.tight_layout()
    return fig


def _persist_weight_draws(df_neural: pd.DataFrame, ranks: np.ndarray,
                          short: pd.DataFrame, out_dir: Path) -> None:
    """Save the full Dirichlet-draw rank data and mark top-10 replacements.

    Writes two CSVs next to the weight-sensitivity figure:

    - ``weight_sensitivity_draws.csv``: one row per (draw, candidate) —
      the rank of every one of the 99 neural candidates under each of the
      1000 weight draws, with an ``in_top_10`` flag;
    - ``weight_sensitivity_top10_challengers.csv``: per-candidate summary
      for everyone who reached the top 10 in at least one draw — baseline
      rank / membership, number and fraction of draws inside the top 10,
      best/median rank, and which baseline Top-10 member was displaced
      most often (and how many times).
    """
    ids = df_neural["gene_id"].astype(str).to_numpy()
    names = (df_neural["gene_name"].astype(str).to_numpy()
             if "gene_name" in df_neural.columns else ids)
    baseline = short["gene_id"].astype(str).tolist()
    base_set = set(baseline)
    base_rank = {gid: i + 1 for i, gid in enumerate(baseline)}
    order = df_neural.sort_values("integrated_score", ascending=False)
    brk = {gid: i + 1 for i, gid in enumerate(order["gene_id"].astype(str))}
    idx_of = {gid: i for i, gid in enumerate(ids)}
    n_draws = ranks.shape[0]

    rows = []
    for d in range(n_draws):
        for i, gid in enumerate(ids):
            r = int(ranks[d, i])
            rows.append({"draw": d + 1, "gene_id": gid,
                         "gene_name": names[i], "rank": r,
                         "in_top_10": r <= 10})
    pd.DataFrame(rows).to_csv(out_dir / "weight_sensitivity_draws.csv",
                              index=False)

    summary = []
    for i, gid in enumerate(ids):
        r = ranks[:, i]
        n_in = int((r <= 10).sum())
        if n_in == 0:
            continue
        displ: dict[str, int] = {}
        for d in range(n_draws):
            if r[d] <= 10:
                for b in baseline:
                    bi = idx_of[b]
                    if ranks[d, bi] > 10:
                        displ[b] = displ.get(b, 0) + 1
        top = max(displ.items(), key=lambda kv: kv[1]) if displ else ("", 0)
        summary.append({
            "gene_id": gid,
            "gene_name": names[i],
            "baseline_rank": brk.get(gid),
            "in_baseline_top10": gid in base_set,
            "n_draws_in_top10": n_in,
            "frac_draws_in_top10": round(n_in / n_draws, 4),
            "best_rank": int(r.min()),
            "median_rank": round(float(np.median(r)), 1),
            "displaced_most_often": top[0],
            "n_times_displaced": top[1],
        })
    out = pd.DataFrame(summary).sort_values(
        ["in_baseline_top10", "n_draws_in_top10"],
        ascending=[False, False])
    out.to_csv(out_dir / "weight_sensitivity_top10_challengers.csv",
               index=False)
    print(f"  wrote weight_sensitivity_draws.csv "
          f"({len(rows)} rows, {n_draws} draws x {len(ids)} candidates)")
    print(f"  wrote weight_sensitivity_top10_challengers.csv "
          f"({len(out)} rows)")


# ---------------------------------------------------------------------------
# PNG-export wrappers (script path). Each calls the corresponding
# `make_*` builder and saves the returned Figure to disk.
# ---------------------------------------------------------------------------

MAIN_FIGURES = {
    "1_score_distributions",
    "2_candidate_summary",
    "3_top10_dual_track",
    "4_evidence_heatmap",
    "5_candidate_funnel",
    "6_evidence_composition",
    "7_stream_ablation",
    "8_top10_radar",
    "9_go_dotplot",
    "10_integrated_vs_composite",
    "11_proof_status_violin",
    "12_weight_sensitivity",
    "13_integrated_vs_neural_filter",
    "14_evidence_heatmap_all_genes",
}


def fig_score_distributions(df: pd.DataFrame, out: Path) -> None:
    f = make_score_distributions(df)
    if f is not None:
        _save(f, out, "1_score_distributions")


def fig_candidate_summary(df: pd.DataFrame, out: Path) -> None:
    f = make_candidate_summary(df)
    if f is not None:
        _save(f, out, "2_candidate_summary")


def fig_top10_dual_track(df_neural: pd.DataFrame | None, out: Path) -> None:
    f = make_top10_dual_track(df_neural)
    if f is not None:
        _save(f, out, "3_top10_dual_track")


def fig_evidence_heatmap(df: pd.DataFrame, out: Path, n: int = 30) -> None:
    f = make_evidence_heatmap(df, n=n)
    if f is not None:
        _save(f, out, "4_evidence_heatmap")


def fig_evidence_heatmap_all_genes(df: pd.DataFrame | None, out: Path) -> None:
    if df is None or df.empty:
        return
    f = make_evidence_heatmap_all_genes(df)
    if f is not None:
        _save(f, out, "14_evidence_heatmap_all_genes")


def fig_candidate_funnel(df: pd.DataFrame, df_neural: pd.DataFrame | None,
                         out: Path) -> None:
    f = make_candidate_funnel(df, df_neural)
    if f is not None:
        _save(f, out, "5_candidate_funnel")


def fig_evidence_composition(df: pd.DataFrame, out: Path, n: int = 15) -> None:
    f = make_evidence_composition(df, n=n)
    if f is not None:
        _save(f, out, "6_evidence_composition")


def fig_stream_ablation(df: pd.DataFrame, out: Path, n: int = 30) -> None:
    f = make_stream_ablation(df, n=n)
    if f is not None:
        _save(f, out, "7_stream_ablation")


def fig_top10_radar(df_neural: pd.DataFrame | None, out: Path) -> None:
    f = make_top10_radar(df_neural)
    if f is not None:
        _save(f, out, "8_top10_radar")


def fig_go_dotplot(df_neural: pd.DataFrame | None, out: Path) -> None:
    f = make_go_dotplot(df_neural)
    if f is not None:
        _save(f, out, "9_go_dotplot")


def fig_integrated_vs_composite(df_neural: pd.DataFrame | None,
                                out: Path) -> None:
    f = make_integrated_vs_composite(df_neural)
    if f is not None:
        _save(f, out, "10_integrated_vs_composite")


def fig_proof_status_violin(df: pd.DataFrame, out: Path) -> None:
    f = make_proof_status_violin(df)
    if f is not None:
        _save(f, out, "11_proof_status_violin")


def fig_weight_sensitivity(df_neural: pd.DataFrame | None, out: Path) -> None:
    f = make_weight_sensitivity(df_neural, out_dir=out)
    if f is not None:
        _save(f, out, "12_weight_sensitivity")


def fig_neural_filter_vs_score(df: pd.DataFrame, out: Path) -> None:
    f = make_neural_filter_vs_score(df)
    if f is not None:
        _save(f, out, "13_integrated_vs_neural_filter")


def main() -> None:
    args = parse_args()
    run_dir: Path = args.run
    out_dir: Path = args.out

    rank_csv = run_dir / "rank.csv"
    neural_csv = run_dir / "rank_neural.csv"

    if not rank_csv.exists():
        raise SystemExit(
            f"rank.csv not found in {run_dir}.\n"
            "Run `bioforge neuraltf run` first, or pass --run <dir>."
        )

    print(f"Loading {rank_csv.name}...")
    df = pd.read_csv(rank_csv)
    print(f"  {len(df)} candidates with {len(df.columns)} columns: {list(df.columns)}")
    df_neural = None
    if neural_csv.exists():
        df_neural = pd.read_csv(neural_csv)
        print(f"  {len(df_neural)} neural-filtered candidates")
    else:
        print("  rank_neural.csv missing — neural-only figures will be skipped")

    print(f"Generating figures into {out_dir}/")
    fig_score_distributions(df, out_dir)
    fig_candidate_summary(df, out_dir)
    fig_top10_dual_track(df_neural, out_dir)
    fig_evidence_heatmap(df, out_dir, n=30)
    fig_candidate_funnel(df, df_neural, out_dir)
    fig_evidence_composition(df, out_dir, n=15)
    fig_stream_ablation(df, out_dir, n=30)
    fig_top10_radar(df_neural, out_dir)
    fig_go_dotplot(df_neural, out_dir)
    fig_integrated_vs_composite(df_neural, out_dir)
    fig_proof_status_violin(df_neural if df_neural is not None else df, out_dir)
    fig_weight_sensitivity(df_neural, out_dir)
    fig_neural_filter_vs_score(df, out_dir)
    fig_evidence_heatmap_all_genes(df_neural if df_neural is not None else df, out_dir)

    stale = [p for p in out_dir.glob("*.png") if p.stem not in MAIN_FIGURES]
    for p in stale:
        p.unlink()
        print(f"  removed stale {p.name}")

    print(f"\nAll figures saved to: {out_dir}")
    print("Generated:")
    for p in sorted(out_dir.glob("*.png")):
        print(f"  {p.name}")


if __name__ == "__main__":
    main()