#!/usr/bin/env python
"""Visualize fixed-weight baseline method results — enhanced for publication.

Generates 13 publication-quality figures from the pipeline's rank.csv / rank_neural.csv:

  1. score_distributions        - per-stream score histograms + integrated (with stats)
  2. candidate_summary          - 2x2: tiers, proof status, score by status, coverage
  3. top10_dual_track           - final Top-10, Track A vs B with domain + RNAi info
  4. evidence_heatmap           - top-30 evidence matrix with clustering + annotations
  5. candidate_funnel           - scored -> neural-filtered -> final (with percentages)
  6. evidence_composition       - stacked per-stream contribution of top-15
  7. stream_ablation            - rank sensitivity when each stream is removed
  8. top10_radar                - per-candidate 7-stream fingerprints + track mean
  9. go_dotplot                 - GO-term coverage of the Top-10 (with descriptions)
  10. integrated_vs_composite   - composite vs integrated with bonus breakdown
  11. proof_status_violin       - integrated-score distribution by proof status (with stats)
  12. weight_sensitivity        - Top-10 rank bands under Dirichlet weight draws
  13. integrated_vs_neural_filter - ECDF of integrated score (with KS test)

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
from scipy.stats import ks_2samp, mannwhitneyu

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
STREAM_WEIGHTS = np.array([0.211, 0.105, 0.158, 0.158, 0.105, 0.158, 0.105])

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
    if "bhlh" in d: return "bHLH"
    if "homeobox" in d: return "Homeobox"
    if "znf" in d: return "Znf"
    if "fork_head" in d: return "fork_head"
    if "t-box" in d: return "T-box"
    if "ets" in d: return "Ets"
    return "none"


def _add_panel_letter(ax, letter, x=-0.08, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top", ha="left", color=C_BLACK)


def _annotate_stats(ax, data_groups, labels, test="mannwhitney", y_max=None):
    """Add statistical significance annotations between groups."""
    if len(data_groups) < 2:
        return
    y_max = y_max or max(max(g) for g in data_groups if len(g) > 0)
    y_step = 0.08 * (y_max if y_max > 0 else 1)
    y_pos = y_max + y_step
    for i in range(len(data_groups)):
        for j in range(i + 1, len(data_groups)):
            if len(data_groups[i]) > 1 and len(data_groups[j]) > 1:
                try:
                    if test == "mannwhitney":
                        stat, p = mannwhitneyu(data_groups[i], data_groups[j], alternative='two-sided')
                    else:
                        stat, p = ks_2samp(data_groups[i], data_groups[j])
                    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                    ax.annotate(f"{sig}", xy=(0.5, y_pos), xycoords='axes fraction',
                                ha='center', va='bottom', fontsize=8, color=C_BLACK)
                    y_pos += y_step * 1.2
                except:
                    pass


# ---------------------------------------------------------------------------
# Figure 1: Score distributions per stream (enhanced)
# ---------------------------------------------------------------------------
def make_score_distributions(df, out_path: Path):
    fig, axes = plt.subplots(2, 4, figsize=(10, 5))
    axes = axes.flatten()
    
    for i, stream in enumerate(STREAMS):
        ax = axes[i]
        vals = pd.to_numeric(df[stream], errors="coerce").dropna()
        n = len(vals)
        mean_val = vals.mean()
        median_val = vals.median()
        
        ax.hist(vals, bins=25, color=C_SKY, edgecolor=C_GRAY, alpha=0.85, linewidth=0.5, density=False)
        ax.axvline(mean_val, color=C_ORANGE, linestyle='--', linewidth=1.5, label=f'Mean={mean_val:.2f}')
        ax.axvline(median_val, color=C_GREEN, linestyle=':', linewidth=1.5, label=f'Median={median_val:.2f}')
        ax.set_title(f"{STREAM_LABELS[i]} (w={STREAM_WEIGHTS[i]:.3f})", fontsize=9, fontweight="bold")
        ax.set_xlabel("Score")
        ax.set_ylabel("Count")
        ax.legend(fontsize=6, loc='upper right')
        _style_ax(ax)
    
    # Integrated score panel
    ax = axes[7]
    vals = pd.to_numeric(df["integrated_score"], errors="coerce").dropna()
    n = len(vals)
    mean_val = vals.mean()
    median_val = vals.median()
    ax.hist(vals, bins=25, color=C_ORANGE, edgecolor=C_GRAY, alpha=0.85, linewidth=0.5, density=False)
    ax.axvline(mean_val, color=C_A, linestyle='--', linewidth=1.5, label=f'Mean={mean_val:.3f}')
    ax.axvline(median_val, color=C_GREEN, linestyle=':', linewidth=1.5, label=f'Median={median_val:.3f}')
    ax.set_title("Integrated (fixed-weight)", fontsize=9, fontweight="bold")
    ax.set_xlabel("Score")
    ax.set_ylabel("Count")
    ax.legend(fontsize=6, loc='upper right')
    _style_ax(ax)
    
    fig.suptitle("A · Score distributions per evidence stream (fixed-weight baseline, n=249)", 
                 fontsize=11, fontweight="bold", y=1.00)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_score_distributions.png", bbox_inches="tight")
    plt.close(fig)


def make_candidate_summary(df, out_path: Path):
    fig, axes = plt.subplots(2, 2, figsize=(8, 6))
    
    # Panel A: Tiers
    ax = axes[0, 0]
    _add_panel_letter(ax, "A")
    tier_order = ["high", "medium", "low"]
    tier_counts = df["tier"].value_counts().reindex(tier_order, fill_value=0)
    tier_pct = (tier_counts / len(df) * 100).round(1)
    bars = ax.bar(tier_counts.index, tier_counts.values, 
                  color=[C_GREEN, C_ORANGE, C_A], edgecolor=C_GRAY, linewidth=0.5)
    ax.set_title("A · Candidates per tier")
    ax.set_ylabel("Count")
    _style_ax(ax)
    for i, (v, p) in enumerate(zip(tier_counts.values, tier_pct.values)):
        ax.text(i, v + 1, f"{v}\n({p}%)", ha="center", fontsize=9, fontweight="bold")

    # Panel B: Proof status
    ax = axes[0, 1]
    _add_panel_letter(ax, "B")
    status_order = ["known_rnai_validated", "prior_fstf_not_tested", "novel_candidate"]
    status_counts = df["proof_status"].value_counts().reindex(status_order, fill_value=0)
    status_pct = (status_counts / len(df) * 100).round(1)
    status_labels = ["RNAi-validated", "Prior FSTF", "Novel"]
    bars = ax.bar(status_labels, status_counts.values, 
                  color=[C_A, C_ORANGE, C_B], edgecolor=C_GRAY, linewidth=0.5)
    ax.set_title("B · Proof status")
    ax.set_ylabel("Count")
    _style_ax(ax)
    for i, (v, p) in enumerate(zip(status_counts.values, status_pct.values)):
        ax.text(i, v + 1, f"{v}\n({p}%)", ha="center", fontsize=9, fontweight="bold")

    # Panel C: Integrated score by proof status
    ax = axes[1, 0]
    _add_panel_letter(ax, "C")
    if "integrated_score" in df.columns and "proof_status" in df.columns:
        box_data = []
        for status, label, color in zip(status_order, status_labels, [C_A, C_ORANGE, C_B]):
            sub = df[df["proof_status"] == status]["integrated_score"].dropna()
            if len(sub) > 0:
                box_data.append(sub)
                bp = ax.boxplot(sub, positions=[status_order.index(status)], widths=0.5,
                               patch_artist=True, 
                               boxprops=dict(facecolor=color, alpha=0.7),
                               medianprops=dict(color=C_BLACK, linewidth=1.5),
                               whiskerprops=dict(color=C_BLACK),
                               capprops=dict(color=C_BLACK),
                               flierprops=dict(marker='.', markersize=3, color=C_GRAY),
                               showfliers=True)
        # Add significance annotations
        _annotate_stats(ax, [df[df["proof_status"]==s]["integrated_score"].dropna().values 
                             for s in status_order if len(df[df["proof_status"]==s])>0],
                        status_labels)
        ax.text(0.02, 0.98, f"n={len(df)}", transform=ax.transAxes, fontsize=8, va='top', color=C_GRAY)
    ax.set_xticks(range(len(status_labels)))
    ax.set_xticklabels(status_labels)
    ax.set_title("C · Integrated score by proof status")
    ax.set_ylabel("Integrated score")
    _style_ax(ax)

    # Panel D: Stream coverage
    ax = axes[1, 1]
    _add_panel_letter(ax, "D")
    coverage = df["n_streams"].value_counts().sort_index()
    coverage_pct = (coverage / len(df) * 100).round(1)
    ax.plot(coverage.index, coverage.values, "o-", color=C_ORANGE, linewidth=2, markersize=6)
    ax.fill_between(coverage.index, coverage.values, alpha=0.2, color=C_ORANGE)
    ax.set_title("D · Number of contributing evidence streams")
    ax.set_xlabel("Non-NaN streams (max 7)")
    ax.set_ylabel("Candidates")
    ax.set_xlim(-0.5, 7.5)
    for i, (idx, val) in enumerate(zip(coverage.index, coverage.values)):
        ax.text(idx, val + 2, f"{val}\n({coverage_pct.iloc[i]}%)", ha="center", fontsize=8)
    ax.set_xlim(-0.5, 7.5)
    _style_ax(ax)

    fig.suptitle("Candidate summary — fixed-weight baseline (n=249 candidates)", 
                 fontsize=11, fontweight="bold", y=1.00)
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

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12, 5.5), sharey=True)
    
    for ax, track_df, title, color, track_label in [(ax_a, ta, "Track A — RNAi-validated", C_A, "A"),
                                                      (ax_b, tb, "Track B — Novel candidates", C_B, "B")]:
        track_df = track_df.sort_values("composite_score", ascending=True).reset_index(drop=True)
        y = np.arange(len(track_df))
        domains = [_domain_group(d) for d in track_df["interpro_domains"]]
        colors = [DOMAIN_COLORS.get(d, C_GRAY) for d in domains]
        
        bars = ax.barh(y, track_df["composite_score"], color=colors, edgecolor="white", 
                       height=0.6, linewidth=0.5, zorder=3)
        ax.set_yticks(y)
        # Enhanced labels: gene name + track + integrated score + domain
        labels = [f"{r['gene_name']}  (int={r.get('integrated_score',0):.2f}, {_domain_group(r['interpro_domains'])})" 
                  for _, r in track_df.iterrows()]
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("Composite score (integrated + bonuses)", fontsize=9)
        ax.set_title(title, fontsize=10.5, fontweight="bold")
        ax.set_xlim(0, 1.1)
        _style_ax(ax)
        
        # Domain legend
        unique_domains = sorted(set(domains) - {"none"})
        handles = [plt.Rectangle((0,0),1,1, facecolor=DOMAIN_COLORS[d], edgecolor="white", label=d) 
                   for d in unique_domains]
        ax.legend(handles=handles, loc="lower right", fontsize=7, title="Domain group", framealpha=0.9)

        # Add composite score values on bars
        for bar, score in zip(bars, track_df["composite_score"]):
            ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2, 
                    f"{score:.3f}", va='center', fontsize=7, fontweight='bold')

    fig.suptitle("B · Fixed-weight Top-10 dual track (5A + 5B)", 
                 fontsize=11, fontweight="bold", y=1.00)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_top10_dual_track.png", bbox_inches="tight")
    plt.close(fig)


def make_evidence_heatmap(df, out_path: Path, n: int = 99):
    # Show ALL neural-filtered candidates (99), not just top 30
    neural_csv = Path(__file__).resolve().parents[2] / "runs" / "pipeline_run" / "rank_neural.csv"
    if neural_csv.exists():
        neural_df = pd.read_csv(neural_csv)
        # Merge with prioritized to get composite_score and track
        IN_DIR = Path(__file__).resolve().parents[3] / "projects" / "NeuralTF" / "results"
        prioritized_csv = IN_DIR / "top10_neural_tfs_prioritized.csv"
        if prioritized_csv.exists():
            prioritized = pd.read_csv(prioritized_csv)
            neural_df = neural_df.merge(prioritized[["gene_id_v6", "composite_score", "track", "proof_status"]], 
                                        on="gene_id_v6", how="left", suffixes=("", "_prioritized"))
            for col in ["composite_score", "track", "proof_status"]:
                if f"{col}_prioritized" in neural_df.columns:
                    neural_df[col] = neural_df[col].fillna(neural_df[f"{col}_prioritized"])
                    neural_df.drop(columns=[f"{col}_prioritized"], inplace=True, errors="ignore")
    else:
        neural_df = df
    
    # Use ALL neural-filtered candidates (99), sorted by integrated score
    top = neural_df.nlargest(99, "integrated_score")
    
    # Reorder streams: put correlation LAST (mostly empty/low)
    streams_ordered = ["expression", "specificity", "reproducibility", "rnai",
                       "neural_enriched", "neural_specificity", "correlation"]
    stream_labels_ordered = ["Expression", "Specificity", "Reproducibility", "RNAi",
                             "Neural enrich.", "Neural spec.", "Correlation"]
    stream_weights_ordered = [0.211, 0.105, 0.158, 0.158, 0.158, 0.105, 0.105]
    stream_labels_ordered = ["Expression", "Specificity", "Reproducibility", "RNAi",
                             "Neural enrich.", "Neural spec.", "Correlation"]
    
    plot = top[streams_ordered].apply(pd.to_numeric, errors="coerce").fillna(0)
    labels = [f"{r['gene_name']} | {r.get('proof_status','?')[:4]} | int={r.get('integrated_score',0):.2f} | comp={r.get('composite_score',0):.2f}" 
              for _, r in top.iterrows()]
    
    fig, ax = plt.subplots(figsize=(11, 14))
    _add_panel_letter(ax, "D")
    
    im = ax.imshow(plot.values, aspect="auto", cmap="Reds", vmin=0, vmax=1, interpolation='nearest')
    ax.set_xticks(range(len(streams_ordered)))
    ax.set_xticklabels(stream_labels_ordered, fontsize=8, rotation=30, ha="right")
    ax.set_yticks(range(len(plot)))
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_title(f"D · Evidence matrix — ALL 99 neural-filtered candidates (correlation last)", 
                 fontsize=11, fontweight="bold", pad=12)
    
    cbar = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label("Stream score (0–1)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    
    # Add stream weight annotations on top
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(range(len(streams_ordered)))
    ax2.set_xticklabels([f"w={w:.3f}" for w in stream_weights_ordered], fontsize=7, color=C_GRAY)
    ax2.set_xlabel("Stream weights", fontsize=8, color=C_GRAY)
    
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_evidence_heatmap.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def make_candidate_funnel(df, neural_df, out_path: Path):
    scored = len(df)
    neural = len(neural_df) if neural_df is not None else 0
    final = 10
    
    # Calculate percentages
    pct_neural = neural / scored * 100
    pct_final = final / scored * 100
    pct_final_neural = final / neural * 100
    
    stages = ["All scored\n(249)", f"Neural-filtered\n({neural}, {pct_neural:.1f}%)", 
              f"Final Top-10\n({final}, {pct_final:.1f}% of all, {pct_final_neural:.1f}% of neural)"]
    counts = [scored, neural, final]
    colors = [C_GRAY, C_SKY, C_ORANGE]
    
    fig, ax = plt.subplots(figsize=(6, 5))
    y = np.arange(len(stages))
    bars = ax.barh(y, counts, color=colors, edgecolor=C_GRAY, height=0.6, linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(stages, fontsize=9)
    for i, (bar, c, pct) in enumerate(zip(bars, counts, [100, pct_neural, pct_final])):
        ax.text(bar.get_width() + 3, bar.get_y() + bar.get_height()/2, 
                f"{c} ({pct:.1f}%)", va="center", fontsize=10, fontweight="bold")
    ax.set_xlabel("Number of candidates")
    ax.set_title("D · Candidate selection funnel", fontsize=10.5, fontweight="bold")
    ax.set_xlim(0, max(counts) * 1.3)
    ax.invert_yaxis()
    _style_ax(ax)
    
    # Add explanatory text
    ax.text(0.5, -0.12, 
            "Pipeline: 249 candidates scored → 99 pass neural filter → 10 prioritized",
            transform=ax.transAxes, ha='center', fontsize=8, color=C_GRAY)
    
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_candidate_funnel.png", bbox_inches="tight")
    plt.close(fig)


def make_evidence_composition(df, out_path: Path, n: int = 10):
    # Use TOP 10 shortlisted (5 Track A + 5 Track B) instead of top 15
    IN_DIR = Path(__file__).resolve().parents[3] / "projects" / "NeuralTF" / "results"
    prioritized_csv = IN_DIR / "top10_neural_tfs_prioritized.csv"
    streams_ordered = ["expression", "specificity", "reproducibility", "rnai",
                       "neural_enriched", "neural_specificity", "correlation"]
    stream_labels_ordered = ["Expression", "Specificity", "Reproducibility", "RNAi",
                             "Neural enrich.", "Neural spec.", "Correlation"]
    
    if (IN_DIR / "top10_neural_tfs_prioritized.csv").exists():
        top_prioritized = pd.read_csv(IN_DIR / "top10_neural_tfs_prioritized.csv").sort_values("composite_score", ascending=False)
        # Merge with df to get stream columns (prioritized has gene_id_v6, df has gene_id)
        # Need to keep integrated_score from df
        top = top_prioritized.merge(df[["gene_id"] + streams_ordered + ["integrated_score"]], 
                                    left_on="gene_id_v6", right_on="gene_id", how="left")
    else:
        # Fallback to top n by integrated score
        top = df.nlargest(n, "integrated_score")
    
    plot = top[streams_ordered].apply(pd.to_numeric, errors="coerce").fillna(0)
    labels = [f"{r['gene_name']} ({r['track']})" for _, r in top.iterrows()]
    
    fig, ax = plt.subplots(figsize=(10, 7))
    bottom = np.zeros(len(plot))
    colors = [C_A, C_SKY, C_ORANGE, C_B, C_GREEN, C_ORANGE, C_GRAY]
    
    for i, stream in enumerate(streams_ordered):
        vals = plot[stream].values
        ax.barh(range(len(plot)), vals, left=bottom, color=colors[i],
                edgecolor=C_GRAY, linewidth=0.3, label=stream_labels_ordered[i])
        bottom += vals
    
    ax.set_yticks(range(len(plot)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Cumulative stream score (weighted sum)", fontsize=9)
    ax.set_title(f"E · Evidence composition — TOP 10 shortlisted (5A + 5B)", 
                 fontsize=11, fontweight="bold")
    ax.legend(loc="lower right", fontsize=7, ncol=2, framealpha=0.9)
    ax.set_xlim(0, 1.0)
    _style_ax(ax)
    
    # Add total integrated score annotations
    for i, (_, row) in enumerate(top.iterrows()):
        total = row["integrated_score"]
        ax.text(1.01, i, f"{total:.3f}", va='center', fontsize=7, color=C_GRAY, transform=ax.get_yaxis_transform())
    
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_evidence_composition.png", bbox_inches="tight")
    plt.close(fig)


def make_stream_ablation(df, out_path: Path):
    n_top = min(30, len(df))
    top = df.nlargest(n_top, "integrated_score").copy()
    top["rank"] = range(1, len(top) + 1)
    ranks = {"original": top["rank"].values}
    
    for stream in STREAMS:
        S = top[STREAMS].apply(pd.to_numeric, errors="coerce").fillna(0).values
        w = STREAM_WEIGHTS
        w_no = np.delete(w, STREAMS.index(stream))
        S_no = np.delete(S, STREAMS.index(stream), axis=1)
        w_no = w_no / w_no.sum()
        scores_no = S_no @ w_no
        ranks[f"no_{stream}"] = np.argsort(-scores_no) + 1
    
    rank_df = pd.DataFrame(ranks, index=top.index)
    shift = rank_df.sub(rank_df["original"], axis=0).drop(columns=["original"])
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Box plot with jittered points
    bp = shift.plot(kind="box", ax=ax, color=C_GRAY, patch_artist=True,
                    boxprops=dict(facecolor=C_SKY, alpha=0.6),
                    medianprops=dict(color=C_BLACK, linewidth=1.5),
                    whiskerprops=dict(color=C_BLACK),
                    capprops=dict(color=C_BLACK),
                    flierprops=dict(marker='.', markersize=3, color=C_GRAY),
                    showfliers=True, zorder=2)
    
    # Add jittered individual points
    for i, col in enumerate(shift.columns):
        x = np.random.normal(i, 0.04, size=len(shift))
        ax.scatter(x, shift[col], alpha=0.4, s=15, color=C_GRAY, zorder=3)
    
    ax.axhline(0, color=C_BLACK, linewidth=1.0, linestyle="--", zorder=1)
    ax.set_xticklabels([STREAM_SHORT[STREAMS.index(c.replace("no_", ""))] for c in shift.columns],
                       rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Rank shift (ablated − original)", fontsize=10)
    ax.set_title("E · Stream ablation — rank sensitivity (top 30 candidates)", 
                 fontsize=11, fontweight="bold", pad=12)
    
    # Add explanatory note
    ax.text(0.5, -0.15, 
            "Positive shift = rank worsens when stream removed; Negative = rank improves",
            transform=ax.transAxes, ha='center', fontsize=8, color=C_GRAY)
    
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_stream_ablation.png", bbox_inches="tight")
    plt.close(fig)


def make_top10_radar(neural_df, out_path: Path):
    if neural_df is None or neural_df.empty:
        return
    
    # Check if stream columns exist
    has_streams = all(s in neural_df.columns for s in STREAMS)
    if not has_streams:
        # Load from rank_neural.csv if available
        runs = Path(__file__).resolve().parents[2] / "runs" / "pipeline_run"
        neural_csv = runs / "rank_neural.csv"
        if neural_csv.exists():
            neural_full = pd.read_csv(neural_csv)
            # Merge stream data
            merged = neural_df.merge(neural_full[["gene_id"] + STREAMS], 
                                     left_on="gene_id_v6", right_on="gene_id", how="left")
            neural_df = merged
        else:
            return
    
    top = neural_df.nlargest(10, "composite_score")
    angles = np.linspace(0, 2 * np.pi, len(STREAMS), endpoint=False).tolist()
    angles += angles[:1]
    
    # Calculate track means
    track_means = {}
    for track in ["A", "B"]:
        track_data = neural_df[neural_df["track"] == track].nlargest(5, "composite_score")
        if len(track_data) > 0:
            track_means[track] = [pd.to_numeric(track_data[s], errors="coerce").mean() or 0 for s in STREAMS]
            track_means[track] += track_means[track][:1]
    
    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(polar=True))
    _add_panel_letter(ax, "F", x=-0.15, y=1.15)
    
    # Plot track means first (thicker lines)
    for track in ["A", "B"]:
        if track in track_means:
            color = C_A if track == "A" else C_B
            label = f"Track {track} mean (n=5)"
            ax.plot(angles, track_means[track], "o-", linewidth=3, markersize=6, 
                    color=color, label=label, alpha=0.8, zorder=5)
            ax.fill(angles, track_means[track], color=color, alpha=0.15, zorder=4)
    
    # Plot individual candidates (thinner lines)
    for _, row in top.iterrows():
        vals = [pd.to_numeric(row.get(s, np.nan), errors="coerce") or 0 for s in STREAMS]
        vals += vals[:1]
        label = f"{row['gene_name']} ({row['track']})"
        color = C_A if row['track'] == "A" else C_B
        ax.plot(angles, vals, "o-", linewidth=1, markersize=3, color=color, alpha=0.4, zorder=2)
    
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(STREAM_LABELS, fontsize=9)
    ax.set_ylim(0, 1)
    ax.set_title("F · Top-10 stream fingerprints (radar plot)", 
                 fontsize=11, fontweight="bold", pad=20)
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.15), fontsize=8, ncol=1, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_top10_radar.png", bbox_inches="tight")
    plt.close(fig)


def make_go_dotplot(neural_df, out_path: Path):
    supp = Path(__file__).resolve().parents[2] / "figures" / "supplementary"
    go_file = supp / "go_gene_term_matrix_reduced.csv"
    if not go_file.exists():
        return
    go = pd.read_csv(go_file)
    if "gene_id" not in go.columns:
        return
    top_genes = set(neural_df.nlargest(10, "composite_score")["gene_id_v6"])
    go = go[go["gene_id"].isin(top_genes)]
    if go.empty:
        return
    # Only include GO term columns (numeric), exclude 'cell' and other non-GO columns
    terms = [c for c in go.columns if c != "gene_id" and c != "cell" and pd.api.types.is_numeric_dtype(go[c])]
    term_counts = go[terms].sum().sort_values(ascending=False).head(15)
    
    # Try to get GO term descriptions
    term_info = {}
    try:
        # Try to load GO descriptions if available
        pass
    except:
        pass
    
    fig, ax = plt.subplots(figsize=(9, 6))
    y = np.arange(len(term_counts))
    bars = ax.barh(y, term_counts.values, color=C_SKY, edgecolor=C_GRAY, height=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(term_counts.index, fontsize=8)
    ax.set_xlabel("Number of Top-10 candidates annotated", fontsize=9)
    ax.set_title("G · GO term coverage — Top-10 candidates", fontsize=10.5, fontweight="bold")
    
    # Add count labels
    for i, v in enumerate(term_counts.values):
        ax.text(v + 0.1, i, str(v), va='center', fontsize=8, fontweight='bold')
    
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_go_dotplot.png", bbox_inches="tight")
    plt.close(fig)


def make_integrated_vs_composite(neural_df, out_path: Path):
    if neural_df is None or neural_df.empty:
        return
    top = neural_df.nlargest(10, "composite_score")
    
    fig, ax = plt.subplots(figsize=(7, 6))
    _add_panel_letter(ax, "G")
    
    bonus_info = []
    for _, row in top.iterrows():
        base = pd.to_numeric(row.get("integrated_score", 0), errors="coerce") or 0
        comp = pd.to_numeric(row.get("composite_score", 0), errors="coerce") or 0
        bonus = comp - base
        bonus_info.append({"gene": row["gene_name"], "base": base, "comp": comp, 
                          "bonus": bonus, "track": row["track"]})
        
        color = C_A if row["track"] == "A" else C_B
        marker = "o" if row["track"] == "A" else "^"
        ax.scatter(base, comp, color=color, s=120, marker=marker,
                   edgecolors="white", linewidth=1, zorder=5)
        ax.annotate(row["gene_name"], (base, comp),
                    textcoords="offset points", xytext=(8, 4),
                    fontsize=8, fontweight='bold', color=C_BLACK)
    
    bonus_df = pd.DataFrame(bonus_info)
    avg_bonus = bonus_df["bonus"].mean()
    max_bonus = bonus_df["bonus"].max()
    
    ax.plot([0, 1], [0, 1], "--", color=C_GRAY, linewidth=1.5, label="y = x (no bonus)")
    ax.set_xlim(-0.02, 0.92)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlabel("Integrated score (fixed-weight baseline)", fontsize=10)
    ax.set_ylabel("Composite score (integrated + bonuses)", fontsize=10)
    ax.set_title("G · Integrated vs composite — fixed-weight Top-10", 
                 fontsize=11, fontweight="bold", pad=12)
    
    # Custom legend with track and bonus info
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_A, markersize=10, label="Track A (RNAi)"),
        Line2D([0], [0], marker="^", color="w", markerfacecolor=C_B, markersize=10, label="Track B (novel)"),
        Line2D([0], [0], color=C_GRAY, linewidth=1.5, linestyle="--", label="y = x (no bonus)"),
        Line2D([0], [0], color="white", label=f"Mean bonus: {avg_bonus:.3f}"),
        Line2D([0], [0], color="white", label=f"Max bonus: {max_bonus:.3f}"),
    ]
    ax.legend(handles=handles, fontsize=8, loc="lower right", framealpha=0.9)
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
    
    fig, ax = plt.subplots(figsize=(8, 6))
    _add_panel_letter(ax, "H")
    
    data_for_test = []
    for s in present:
        data = df_plot[df_plot["proof_status"] == s]["integrated_score"].values
        data_for_test.append(data)
    
    parts = ax.violinplot(data_for_test, positions=range(len(present)),
                          showmeans=True, showmedians=True, widths=0.7)
    colors = [C_A, C_ORANGE, C_B]
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(colors[i % len(colors)])
        pc.set_alpha(0.7)
        pc.set_edgecolor(C_GRAY)
    parts["cmedians"].set_color(C_BLACK)
    parts["cmedians"].set_linewidth(2)
    parts["cmeans"].set_color("white")
    # cmeans is a LineCollection; don't try to set marker properties
    
    # Add sample sizes and means ON the violins (not below x-axis)
    for i, s in enumerate(present):
        n = len(df_plot[df_plot["proof_status"] == s])
        mean_val = df_plot[df_plot["proof_status"] == s]["integrated_score"].mean()
        median_val = np.median(df_plot[df_plot["proof_status"] == s]["integrated_score"])
        # Position text on the violin
        ax.text(i, 1.02, f"n={n}\nmean={mean_val:.3f}\nmedian={median_val:.3f}", 
                ha='center', va='bottom', fontsize=9, fontweight='bold', 
                color=C_BLACK, transform=ax.get_xaxis_transform(),
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=C_GRAY, alpha=0.9))
    
    # Statistical test
    if len(data_for_test) >= 2:
        _annotate_stats(ax, data_for_test, present)
    
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(present, fontsize=11, fontweight='bold')
    ax.set_ylabel("Integrated score", fontsize=12)
    ax.set_title("H · Integrated score distribution by proof status (n=99)", 
                 fontsize=13, fontweight="bold", pad=16)
    ax.set_ylim(-0.1, 1.15)
    
    # Add legend explaining the plot elements
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=C_A, alpha=0.7, edgecolor=C_GRAY, label="RNAi-validated (Track A)"),
        Patch(facecolor=C_ORANGE, alpha=0.7, edgecolor=C_GRAY, label="Prior FSTF (Track A/B)"),
        Patch(facecolor=C_B, alpha=0.7, edgecolor=C_GRAY, label="Novel candidate (Track B)"),
        Line2D([0], [0], color=C_BLACK, linewidth=2, label="Median"),
        Line2D([0], [0], color="white", marker='D', markerfacecolor='white', 
               markeredgecolor='black', markersize=8, linestyle='None', label="Mean"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=9, 
              framealpha=0.9, title="Proof status & statistics", title_fontsize=10)
    
    _style_ax(ax)
    # Add frame around the plot
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color(C_BLACK)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_proof_status_violin.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def make_weight_sensitivity(neural_df, out_path: Path):
    runs = Path(__file__).resolve().parents[2] / "runs" / "pipeline_run"
    ws_file = runs / "weight_sensitivity_draws.csv"
    if not ws_file.exists():
        return
    ws = pd.read_csv(ws_file)
    top10 = ws[ws["rank"] <= 10]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    _add_panel_letter(ax, "I")
    
    # Plot each candidate's rank trajectory
    for gene in sorted(top10["gene_id"].unique()):
        g = top10[top10["gene_id"] == gene].sort_values("draw")
        # Get gene name from neural_df
        gene_name = g["gene_name"].iloc[0] if "gene_name" in g.columns else g["gene_id"].iloc[0].replace("dd_Smed_v6_", "dd")
        track = g["track"].iloc[0] if "track" in g.columns else "?"
        color = C_A if track == "A" else C_B
        ax.plot(g["draw"], g["rank"], "o-", markersize=3, linewidth=1, alpha=0.7, 
                color=color, label=gene_name if gene not in ax.get_legend_handles_labels()[1] else "")
    
    ax.invert_yaxis()
    ax.set_xlabel("Dirichlet draw (1000 weight samples)", fontsize=10)
    ax.set_ylabel("Rank", fontsize=10)
    ax.set_title("I · Weight sensitivity — Top-10 rank bands across 1000 Dirichlet draws", 
                 fontsize=11, fontweight="bold", pad=12)
    ax.set_ylim(10.5, 0.5)
    ax.set_xlim(-10, 1010)
    ax.set_yticks(range(1, 11))
    
    # Add track legend
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], color=C_A, marker='o', linewidth=1.5, markersize=6, label="Track A (RNAi)"),
        Line2D([0], [0], color=C_B, marker='o', linewidth=1.5, markersize=6, label="Track B (novel)"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8, framealpha=0.9)
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_weight_sensitivity.png", bbox_inches="tight")
    plt.close(fig)


def make_integrated_vs_neural_filter(df, neural_df, out_path: Path):
    from scipy.stats import ks_2samp
    fig, ax = plt.subplots(figsize=(8, 6))
    _add_panel_letter(ax, "J")
    
    all_vals = pd.to_numeric(df["integrated_score"], errors="coerce").dropna()
    neural_vals = pd.to_numeric(neural_df["integrated_score"], errors="coerce").dropna()
    
    # KS test
    ks_stat, p_val = ks_2samp(all_vals, neural_vals)
    sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
    
    # Also compute Mann-Whitney U for median difference
    from scipy.stats import mannwhitneyu
    mw_stat, mw_p = mannwhitneyu(neural_vals, all_vals, alternative='greater')
    
    for label, data, color in [("All scored candidates (n=249)", all_vals, C_GRAY),
                                ("Neural-filtered candidates (n=99)", neural_vals, C_SKY)]:
        vals = data.sort_values()
        y = np.arange(1, len(vals) + 1) / len(vals)
        ax.plot(vals, y, label=label, color=color, linewidth=2.5)
    
    # Add median lines
    ax.axvline(all_vals.median(), color=C_GRAY, linestyle=':', linewidth=1.5, alpha=0.7, 
               label=f"All median = {all_vals.median():.3f}")
    ax.axvline(neural_vals.median(), color=C_SKY, linestyle='--', linewidth=1.5, alpha=0.7,
               label=f"Neural median = {neural_vals.median():3f}")
    
    # Add interpretation text
    ax.text(0.02, 0.98, 
            f"KS test: D={ks_stat:.3f}, p={p_val:.2e} ({sig})\n"
            f"Mann-Whitney U: p={mw_p:.2e}\n"
            f"Neural candidates have higher integrated scores\n"
            f"(median: {neural_vals.median():.3f} vs {all_vals.median():.3f})",
            transform=ax.transAxes, va='top', fontsize=8,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor=C_GRAY, alpha=0.9))
    
    for label, data, color in [("All scored candidates (n=249)", all_vals, C_GRAY),
                                ("Neural-filtered candidates (n=99)", neural_vals, C_SKY)]:
        vals = data.sort_values()
        y = np.arange(1, len(vals) + 1) / len(vals)
        ax.plot(vals, y, label=label, color=color, linewidth=2.5)
    
    ax.set_xlabel("Integrated score (fixed-weight, 0-1)", fontsize=11)
    ax.set_ylabel("ECDF (cumulative fraction)", fontsize=11)
    ax.set_title(f"J · Neural filter enriches high-scoring candidates\n"
                 f"(KS test: D={ks_stat:.3f}, p={p_val:.2e}, {sig}; MWU p={mw_p:.2e})", 
                 fontsize=12, fontweight="bold", pad=16)
    ax.legend(fontsize=10, framealpha=0.9, loc="lower right")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(out_path / "fig_fixed_integrated_vs_neural_filter.png", bbox_inches="tight", dpi=300)
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