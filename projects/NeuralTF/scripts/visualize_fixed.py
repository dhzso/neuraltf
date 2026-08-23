#!/usr/bin/env python
"""Visualize fixed-weight baseline method results.

Generates 13 publication-quality figures from the pipeline's rank.csv / rank_neural.csv:

  1. score_distributions        - per-stream score histograms + integrated
  2. candidate_summary          - 2x2: tiers, proof status, score by status, coverage
  3. top10_dual_track           - final Top-10, Track A (RNAi-validated) vs Track B (novel)
  4. evidence_heatmap           - top-30 evidence matrix (core figure)
  5. candidate_funnel           - scored -> neural-filtered -> final candidates
  6. evidence_composition       - stacked per-stream contribution of top-15
  7. stream_ablation            - rank sensitivity when each stream is removed
  8. top10_radar                - per-candidate 7-stream fingerprints
  9. go_dotplot                 - GO-term coverage of the Top-10
  10. integrated_vs_composite   - composite (integrated + bonuses) vs integrated
  11. proof_status_violin       - integrated-score distribution by proof status
  12. weight_sensitivity        - Top-10 rank bands under Dirichlet weight draws
  13. integrated_vs_neural_filter - ECDF of integrated score for all scored candidates

Outputs: projects/NeuralTF/figures/fig_fixed_*.png

Usage:
    python projects/NeuralTF/scripts/visualize_fixed.py
    python projects/NeuralTF/scripts/visualize_fixed.py --run projects/NeuralTF/runs/pipeline_run --out projects/NeuralTF/figures
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
IN_DIR = REPO_ROOT / "projects" / "NeuralTF" / "results"


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

STREAMS = ["expression", "specificity", "reproducibility", "rnai",
           "correlation", "neural_enriched", "neural_specificity"]
STREAM_LABELS = ["Expression", "Specificity", "Reproducibility", "RNAi",
                 "Correlation", "Neural enrich.", "Neural spec."]
STREAM_SHORT = ["Expr.", "Spec.", "Reprod.", "RNAi", "Corr.", "N.enr.", "N.spec."]

DOMAIN_COLORS = {
    "bHLH": C_B,
    "Homeobox": C_GREEN,
    "Znf": "#CC79A7",
    "fork_head": C_ORANGE,
    "T-box": C_A,
    "Ets": C_SKY,
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


def _domain_group(domains: str) -> str:
    if not isinstance(domains, str) or not domains.strip():
        return "none"
    d = domains.lower()
    # Map to exact DOMAIN_COLORS keys
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


# ---------------------------------------------------------------------------
# Figure builders
# ---------------------------------------------------------------------------
def make_score_distributions(df, out_path: Path):
    fig, axes = plt.subplots(2, 4, figsize=(10, 5))
    axes = axes.flatten()
    for i, stream in enumerate(STREAMS):
        ax = axes[i]
        vals = pd.to_numeric(df[stream], errors="coerce").dropna()
        ax.hist(vals, bins=25, color=C_SKY, edgecolor=C_GRAY, alpha=0.85, linewidth=0.5)
        ax.set_title(STREAM_LABELS[i], fontsize=9.5, fontweight="bold")
        ax.set_xlabel("Score")
        ax.set_ylabel("Count")
        _style_ax(ax)
    # Integrated score in last panel
    ax = axes[7]
    vals = pd.to_numeric(df["integrated_score"], errors="coerce").dropna()
    ax.hist(vals, bins=25, color=C_ORANGE, edgecolor=C_GRAY, alpha=0.85, linewidth=0.5)
    ax.set_title("Integrated", fontsize=9.5, fontweight="bold")
    ax.set_xlabel("Score")
    ax.set_ylabel("Count")
    _style_ax(ax)
    fig.suptitle("Score distributions per evidence stream", fontsize=11, fontweight="bold", y=1.00)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_score_distributions.png", bbox_inches="tight")
    plt.close(fig)


def make_candidate_summary(df, out_path: Path):
    fig, axes = plt.subplots(2, 2, figsize=(8, 6))
    ax = axes[0, 0]
    tier_order = ["high", "medium", "low"]
    tier_counts = df["tier"].value_counts().reindex(tier_order, fill_value=0)
    ax.bar(tier_counts.index, tier_counts.values, color=[C_GREEN, C_ORANGE, C_A], edgecolor=C_GRAY, linewidth=0.5)
    ax.set_title("A · Candidates per tier")
    ax.set_ylabel("Count")
    _style_ax(ax)
    for i, v in enumerate(tier_counts.values):
        ax.text(i, v + 0.5, str(v), ha="center", fontsize=9)

    ax = axes[0, 1]
    status_order = ["known_rnai_validated", "prior_fstf_not_tested", "novel_candidate"]
    status_counts = df["proof_status"].value_counts().reindex(status_order, fill_value=0)
    status_labels = ["RNAi-validated", "Prior FSTF", "Novel"]
    ax.bar(status_labels, status_counts.values, color=[C_A, C_ORANGE, C_B], edgecolor=C_GRAY, linewidth=0.5)
    ax.set_title("B · Proof status")
    ax.set_ylabel("Count")
    _style_ax(ax)
    for i, v in enumerate(status_counts.values):
        ax.text(i, v + 0.5, str(v), ha="center", fontsize=9)

    ax = axes[1, 0]
    if "integrated_score" in df.columns and "proof_status" in df.columns:
        for status, label, color in zip(status_order, status_labels, [C_A, C_ORANGE, C_B]):
            sub = df[df["proof_status"] == status]["integrated_score"]
            if len(sub) > 0:
                ax.boxplot(sub, positions=[status_order.index(status)], widths=0.5,
                           patch_artist=True, boxprops=dict(facecolor=color, alpha=0.7),
                           medianprops=dict(color=C_BLACK), showfliers=False)
    ax.set_xticks(range(len(status_labels)))
    ax.set_xticklabels(status_labels)
    ax.set_title("C · Integrated score by proof status")
    ax.set_ylabel("Integrated score")
    _style_ax(ax)

    ax = axes[1, 1]
    coverage = df["n_streams"].value_counts().sort_index()
    ax.plot(coverage.index, coverage.values, "o-", color=C_ORANGE, linewidth=1.5, markersize=5)
    ax.set_title("D · Streams contributing per candidate")
    ax.set_xlabel("Number of non-NaN streams")
    ax.set_ylabel("Candidates")
    ax.set_xlim(0, 7)
    _style_ax(ax)

    fig.suptitle("Candidate summary — fixed-weight baseline", fontsize=11, fontweight="bold", y=1.00)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_candidate_summary.png", bbox_inches="tight")
    plt.close(fig)


def make_top10_dual_track(neural_df, out_path: Path):
    if neural_df is None or neural_df.empty:
        return
    ta = neural_df[neural_df["proof_status"] == "known_rnai_validated"].nlargest(5, "composite_score")
    tb = neural_df[neural_df["proof_status"] == "novel_candidate"].nlargest(5, "composite_score")
    if len(ta) < 5:
        tb = neural_df[neural_df["proof_status"] != "known_rnai_validated"].nlargest(10 - len(ta), "composite_score")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    for ax, track_df, title, color in [(ax_a, ta, "Track A — RNAi-validated", C_A),
                                        (ax_b, tb, "Track B — Novel candidates", C_B)]:
        track_df = track_df.sort_values("composite_score", ascending=True).reset_index(drop=True)
        y = np.arange(len(track_df))
        domains = [_domain_group(d) for d in track_df["interpro_domains"]]
        colors = [DOMAIN_COLORS.get(d, C_GRAY) for d in domains]
        ax.barh(y, track_df["composite_score"], color=colors, edgecolor=C_GRAY, height=0.6, linewidth=0.5)
        ax.set_yticks(y)
        ax.set_yticklabels(track_df["gene_name"], fontsize=9)
        ax.set_xlabel("Composite score")
        ax.set_title(title, fontsize=10.5, fontweight="bold")
        ax.set_xlim(0, 1.05)
        _style_ax(ax)
        # Legend for domains
        handles = [plt.Rectangle((0, 0), 1, 1, facecolor=DOMAIN_COLORS[d], edgecolor=C_GRAY, label=d)
                   for d in sorted(set(domains)) if d != "none"]
        ax.legend(handles=handles, loc="lower right", fontsize=7, title="Domain group")

    fig.suptitle("Fixed-weight Top-10 dual track", fontsize=11, fontweight="bold", y=1.00)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_top10_dual_track.png", bbox_inches="tight")
    plt.close(fig)


def make_evidence_heatmap(df, out_path: Path, n: int = 30):
    top = df.nlargest(n, "integrated_score")
    plot = top[STREAMS].apply(pd.to_numeric, errors="coerce").fillna(0)
    labels = [f"{r['gene_name']} ({r['proof_status'][:3]})" for _, r in top.iterrows()]

    fig, ax = plt.subplots(figsize=(8, 10))
    im = ax.imshow(plot.values, aspect="auto", cmap="Reds", vmin=0, vmax=1)
    ax.set_xticks(range(len(STREAMS)))
    ax.set_xticklabels(STREAM_SHORT, fontsize=8.5, rotation=30, ha="right")
    ax.set_yticks(range(len(plot)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_title(f"Evidence matrix — top {n} by integrated score", fontsize=10.5, fontweight="bold", pad=12)
    plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Stream score (0–1)")
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_evidence_heatmap.png", bbox_inches="tight")
    plt.close(fig)


def make_candidate_funnel(df, neural_df, out_path: Path):
    scored = len(df)
    neural = len(neural_df) if neural_df is not None else 0
    final = 10
    stages = ["All scored", "Neural-filtered", "Final Top-10"]
    counts = [scored, neural, final]

    fig, ax = plt.subplots(figsize=(5, 4))
    y = np.arange(len(stages))
    ax.barh(y, counts, color=[C_GRAY, C_SKY, C_ORANGE], edgecolor=C_GRAY, height=0.6, linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(stages, fontsize=9.5)
    for i, c in enumerate(counts):
        ax.text(c + 2, i, str(c), va="center", fontsize=10, fontweight="bold")
    ax.set_xlabel("Candidates")
    ax.set_title("Candidate funnel", fontsize=10.5, fontweight="bold")
    ax.set_xlim(0, max(counts) * 1.2)
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_candidate_funnel.png", bbox_inches="tight")
    plt.close(fig)


def make_evidence_composition(df, out_path: Path, n: int = 15):
    top = df.nlargest(n, "integrated_score")
    plot = top[STREAMS].apply(pd.to_numeric, errors="coerce").fillna(0)
    labels = [r["gene_name"] for _, r in top.iterrows()]

    fig, ax = plt.subplots(figsize=(8, 8))
    bottom = np.zeros(len(plot))
    colors = [C_A, C_SKY, C_ORANGE, C_B, C_GREEN, C_ORANGE, C_GRAY]
    for i, stream in enumerate(STREAMS):
        ax.barh(range(len(plot)), plot[stream].values, left=bottom, color=colors[i],
                edgecolor=C_GRAY, linewidth=0.3, label=STREAM_LABELS[i])
        bottom += plot[stream].values
    ax.set_yticks(range(len(plot)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Cumulative stream score")
    ax.set_title(f"Evidence composition — top {n}", fontsize=10.5, fontweight="bold")
    ax.legend(loc="lower right", fontsize=7, ncol=2)
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_evidence_composition.png", bbox_inches="tight")
    plt.close(fig)


def make_stream_ablation(df, out_path: Path):
    from scipy.stats import spearmanr
    n_top = min(30, len(df))
    top = df.nlargest(n_top, "integrated_score").copy()
    top["rank"] = range(1, len(top) + 1)
    ranks = {"original": top["rank"].values}
    for stream in STREAMS:
        # Recompute integrated score without this stream
        S = top[STREAMS].apply(pd.to_numeric, errors="coerce").fillna(0).values
        w = np.array([0.211, 0.105, 0.158, 0.158, 0.105, 0.158, 0.105])
        # Remove stream weight and renormalize remaining
        w_no = np.delete(w, STREAMS.index(stream))
        S_no = np.delete(S, STREAMS.index(stream), axis=1)
        w_no = w_no / w_no.sum()
        scores_no = S_no @ w_no
        ranks[f"no_{stream}"] = np.argsort(-scores_no) + 1
    rank_df = pd.DataFrame(ranks, index=top.index)
    shift = rank_df.sub(rank_df["original"], axis=0).drop(columns=["original"])

    fig, ax = plt.subplots(figsize=(8, 6))
    shift.plot(kind="box", ax=ax, color=C_GRAY, patch_artist=True,
               boxprops=dict(facecolor=C_SKY, alpha=0.5),
               medianprops=dict(color=C_BLACK), flierprops=dict(marker=".", markersize=3))
    ax.axhline(0, color=C_BLACK, linewidth=0.8, linestyle="--")
    ax.set_xticklabels([STREAM_SHORT[STREAMS.index(c.replace("no_", ""))] for c in shift.columns],
                       rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Rank shift (ablated − original)")
    ax.set_title("Stream ablation — rank sensitivity (top 30)", fontsize=10.5, fontweight="bold")
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_stream_ablation.png", bbox_inches="tight")
    plt.close(fig)


def make_top10_radar(neural_df, out_path: Path):
    if neural_df is None or neural_df.empty:
        return
    top = neural_df.nlargest(10, "composite_score")
    angles = np.linspace(0, 2 * np.pi, len(STREAMS), endpoint=False).tolist()
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))
    for _, row in top.iterrows():
        vals = [pd.to_numeric(row.get(s, np.nan), errors="coerce") or 0 for s in STREAMS]
        vals += vals[:1]
        label = f"{row['gene_name']} ({row['track']})"
        ax.plot(angles, vals, "o-", linewidth=1.5, markersize=4, label=label)
        ax.fill(angles, vals, alpha=0.1)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(STREAM_SHORT, fontsize=8)
    ax.set_ylim(0, 1)
    ax.set_title("Top-10 stream fingerprints", fontsize=10.5, fontweight="bold", pad=15)
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_top10_radar.png", bbox_inches="tight")
    plt.close(fig)


def make_go_dotplot(neural_df, out_path: Path):
    # Requires PlanMine GO annotations - load from supplementary figures
    supp = Path(__file__).resolve().parents[2] / "figures" / "supplementary"
    go_file = supp / "go_gene_term_matrix_reduced.csv"
    if not go_file.exists():
        return
    go = pd.read_csv(go_file)
    if "gene_id" not in go.columns:
        return
    top_genes = set(neural_df.nlargest(10, "composite_score")["gene_id"])
    go = go[go["gene_id"].isin(top_genes)]
    if go.empty:
        return
    terms = [c for c in go.columns if c != "gene_id"]
    term_counts = go[terms].sum().sort_values(ascending=False).head(15)

    fig, ax = plt.subplots(figsize=(8, 5))
    y = np.arange(len(term_counts))
    ax.barh(y, term_counts.values, color=C_SKY, edgecolor=C_GRAY, height=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(term_counts.index, fontsize=8)
    ax.set_xlabel("Number of Top-10 candidates annotated")
    ax.set_title("GO term coverage — Top-10 candidates", fontsize=10.5, fontweight="bold")
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_go_dotplot.png", bbox_inches="tight")
    plt.close(fig)


def make_integrated_vs_composite(neural_df, out_path: Path):
    if neural_df is None or neural_df.empty:
        return
    top = neural_df.nlargest(10, "composite_score")
    fig, ax = plt.subplots(figsize=(6, 5))
    for _, row in top.iterrows():
        base = pd.to_numeric(row.get("integrated_score", 0), errors="coerce") or 0
        comp = pd.to_numeric(row.get("composite_score", 0), errors="coerce") or 0
        color = C_A if row["track"] == "A" else C_B
        marker = "o" if row["track"] == "A" else "^"
        ax.scatter(base, comp, color=color, s=80, marker=marker,
                   edgecolors="white", linewidth=0.8, zorder=3,
                   label=row["track"] if row["track"] not in ax.get_legend_handles_labels()[1] else "")
        ax.annotate(row["gene_name"], (base, comp),
                    textcoords="offset points", xytext=(5, 3),
                    fontsize=7, color=C_BLACK)
    ax.plot([0, 1], [0, 1], "--", color=C_GRAY, linewidth=1, label="y = x (no bonus)")
    ax.set_xlim(0, 0.9)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Integrated score (fixed-weight baseline)")
    ax.set_ylabel("Composite score (integrated + bonuses)")
    ax.set_title("Integrated vs composite — fixed-weight Top-10", fontsize=10.5, fontweight="bold")
    ax.legend(fontsize=8)
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_integrated_vs_composite.png", bbox_inches="tight")
    plt.close(fig)


def make_proof_status_violin(neural_df, out_path: Path):
    if neural_df is None or neural_df.empty:
        return
    df_plot = neural_df.copy()
    df_plot["proof_status"] = df_plot["proof_status"].replace({
        "known_rnai_validated": "RNAi-validated",
        "prior_fstf_not_tested": "Prior FSTF",
        "novel_candidate": "Novel"
    })
    order = ["RNAi-validated", "Prior FSTF", "Novel"]
    present = [s for s in order if s in df_plot["proof_status"].values]

    fig, ax = plt.subplots(figsize=(6, 4))
    parts = ax.violinplot([df_plot[df_plot["proof_status"] == s]["integrated_score"].values
                           for s in present], positions=range(len(present)),
                          showmeans=True, showmedians=True, widths=0.6)
    colors = [C_A, C_ORANGE, C_B]
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(colors[i % len(colors)])
        pc.set_alpha(0.6)
        pc.set_edgecolor(C_GRAY)
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(present, fontsize=9)
    ax.set_ylabel("Integrated score")
    ax.set_title("Integrated score by proof status", fontsize=10.5, fontweight="bold")
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_proof_status_violin.png", bbox_inches="tight")
    plt.close(fig)


def make_weight_sensitivity(neural_df, out_path: Path):
    # Uses dirichlet_method_comparison's draw data if available
    runs = Path(__file__).resolve().parents[2] / "runs" / "pipeline_run"
    ws_file = runs / "weight_sensitivity_draws.csv"
    if not ws_file.exists():
        return
    ws = pd.read_csv(ws_file)
    top10 = ws[ws["rank"] <= 10]

    fig, ax = plt.subplots(figsize=(8, 5))
    for gene in top10["gene_id"].unique():
        g = top10[top10["gene_id"] == gene]
        ax.plot(g["draw"], g["rank"], "o-", markersize=3, linewidth=0.7, alpha=0.6)
    ax.invert_yaxis()
    ax.set_xlabel("Dirichlet draw")
    ax.set_ylabel("Rank")
    ax.set_title("Weight sensitivity — Top-10 rank bands", fontsize=10.5, fontweight="bold")
    ax.set_ylim(10.5, 0.5)
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_weight_sensitivity.png", bbox_inches="tight")
    plt.close(fig)


def make_integrated_vs_neural_filter(df, neural_df, out_path: Path):
    from scipy.stats import ks_2samp
    fig, ax = plt.subplots(figsize=(6, 4))
    for label, data, color in [("All scored", df["integrated_score"], C_GRAY),
                                ("Neural-filtered", neural_df["integrated_score"], C_SKY)]:
        vals = pd.to_numeric(data, errors="coerce").dropna().sort_values()
        y = np.arange(1, len(vals) + 1) / len(vals)
        ax.plot(vals, y, label=label, color=color, linewidth=2)
    ax.set_xlabel("Integrated score")
    ax.set_ylabel("ECDF")
    ax.set_title("Integrated score distribution — neural filter", fontsize=10.5, fontweight="bold")
    ax.legend(fontsize=8)
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_integrated_vs_neural_filter.png", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    args = parse_args()
    run_dir = args.run
    out_dir = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    rank_csv = run_dir / "rank.csv"
    neural_csv = run_dir / "rank_neural.csv"

    if not rank_csv.exists():
        print(f"[error] {rank_csv} not found")
        return 1

    df = pd.read_csv(rank_csv)
    neural_df = pd.read_csv(neural_csv) if neural_csv.exists() else None

    print(f"Generating fixed-weight figures from {run_dir} -> {out_dir}")
    print(f"  All candidates: {len(df)}, Neural: {len(neural_df) if neural_df is not None else 0}")

    make_score_distributions(df, out_dir)
    print("  wrote fig_fixed_score_distributions.png")
    make_candidate_summary(df, out_dir)
    print("  wrote fig_fixed_candidate_summary.png")
    # Load fixed prioritized CSV for top-10 track data (has composite_score)
    fixed_top10_path = IN_DIR / "top10_neural_tfs_prioritized.csv"
    if fixed_top10_path.exists():
        fixed_top10 = pd.read_csv(fixed_top10_path)
        make_top10_dual_track(fixed_top10, out_dir)
    else:
        make_top10_dual_track(neural_df, out_dir)
    print("  wrote fig_fixed_top10_dual_track.png")
    make_evidence_heatmap(df, out_dir)
    print("  wrote fig_fixed_evidence_heatmap.png")
    make_candidate_funnel(df, neural_df, out_dir)
    print("  wrote fig_fixed_candidate_funnel.png")
    make_evidence_composition(df, out_dir)
    print("  wrote fig_fixed_evidence_composition.png")
    make_stream_ablation(df, out_dir)
    print("  wrote fig_fixed_stream_ablation.png")
    # Use fixed prioritized CSV for top-10 visualizations
    fixed_top10 = pd.read_csv(IN_DIR / "top10_neural_tfs_prioritized.csv") if (IN_DIR / "top10_neural_tfs_prioritized.csv").exists() else neural_df
    make_top10_radar(fixed_top10, out_dir)
    print("  wrote fig_fixed_top10_radar.png")
    make_go_dotplot(fixed_top10, out_dir)
    print("  wrote fig_fixed_go_dotplot.png")
    make_integrated_vs_composite(fixed_top10, out_dir)
    print("  wrote fig_fixed_integrated_vs_composite.png")
    make_proof_status_violin(neural_df, out_dir)
    print("  wrote fig_fixed_proof_status_violin.png")
    make_weight_sensitivity(fixed_top10, out_dir)
    print("  wrote fig_fixed_weight_sensitivity.png")
    make_integrated_vs_neural_filter(df, neural_df, out_dir)
    print("  wrote fig_fixed_integrated_vs_neural_filter.png")

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())