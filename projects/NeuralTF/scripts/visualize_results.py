#!/usr/bin/env python
"""Generate NeuralTF pipeline visualization figures from a run's outputs.

This script reads only what the pipeline actually writes:
  - rank.csv (per-TF integrated scores + 8 evidence streams + proof_status + tier)
  - rank_neural.csv (neural-filtered subset)
  - pipeline_results.json (top 50 candidates with tier)

It produces PNG figures under projects/NeuralTF/figures/ (by default) or
a user-specified output directory. The script is portable — no hardcoded
absolute paths, no hardcoded gene-ID lists; everything is derived from the
input CSVs.

Usage:
    python projects/NeuralTF/scripts/visualize_results.py
    python projects/NeuralTF/scripts/visualize_results.py \\
        --run projects/NeuralTF/runs/pipeline_run \\
        --out projects/NeuralTF/figures
"""
from __future__ import annotations

import argparse
import json
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


# Style ---------------------------------------------------------------------
sns.set_style("whitegrid")
plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})
TIER_COLORS = {"high": "#d62728", "medium": "#ff7f0e", "low": "#2ca02c"}
PROOF_COLORS = {
    "known_rnai_validated": "#1f77b4",
    "prior_fstf_not_tested": "#ff7f0e",
    "novel_candidate": "#2ca02c",
}
STREAM_COLORS = {
    "expression": "#1f77b4",
    "specificity": "#ff7f0e",
    "reproducibility": "#2ca02c",
    "rnai": "#d62728",
    "correlation": "#9467bd",
    "function": "#8c564b",
    "neural_enriched": "#e377c2",
    "neural_specificity": "#7f7f7f",
}
SCORE_STREAMS = list(STREAM_COLORS.keys())


def _tier_color(tier_name: str) -> str:
    return TIER_COLORS.get(str(tier_name).lower(), "#7f7f7f")


def _save(fig, out_dir: Path, name: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"{name}.png"
    fig.tight_layout()
    fig.savefig(p)
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


# ---------------------------------------------------------------------------
# Figure builders — return matplotlib Figure so the same function can be
# invoked by both this PNG-export script and the Streamlit UI (which can
# st.pyplot() a returned Figure directly without re-saving to disk).
# ---------------------------------------------------------------------------

def make_score_distributions(df: pd.DataFrame) -> plt.Figure:
    streams_present = [s for s in SCORE_STREAMS if s in df.columns]
    streams_present += ["integrated_score"]
    n = len(streams_present)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3))
    axes = np.atleast_1d(axes).ravel()
    for ax in axes[n:]:
        ax.set_visible(False)
    for ax, s in zip(axes, streams_present):
        vals = df[s].dropna()
        if vals.empty:
            ax.set_visible(False)
            continue
        ax.hist(vals, bins=20, edgecolor="white", alpha=0.8)
        ax.set_title(s.replace("_", " ").title())
        ax.set_xlabel("Score")
        ax.set_ylabel("Count")
    fig.suptitle(
        f"Evidence Stream Score Distributions ({len(df)} candidates)",
        fontsize=14,
    )
    return fig


def make_tier_distribution(df: pd.DataFrame) -> plt.Figure | None:
    if "tier" not in df.columns:
        return None
    counts = df["tier"].value_counts()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    pal = [_tier_color(t) for t in counts.index]
    axes[0].pie(
        counts.values, labels=counts.index, colors=pal, autopct="%1.0f%%"
    )
    axes[0].set_title("Tier Distribution")
    axes[1].bar(counts.index, counts.values, color=pal, edgecolor="white")
    axes[1].set_ylabel("Number of Candidates")
    axes[1].set_title("Tier Counts")
    return fig


def make_proof_status(df: pd.DataFrame) -> plt.Figure | None:
    if "proof_status" not in df.columns:
        return None
    counts = df["proof_status"].value_counts()
    fig, ax = plt.subplots(figsize=(7, max(3, 0.5 * len(counts))))
    pal = [PROOF_COLORS.get(s, "#999999") for s in counts.index]
    ax.barh(counts.index, counts.values, color=pal, edgecolor="white")
    for i, v in enumerate(counts.values):
        ax.text(v + 0.5, i, str(v), va="center")
    ax.set_xlabel("Number of Candidates")
    ax.set_title("Proof Status Distribution")
    return fig


def make_top_candidates(df: pd.DataFrame, n: int = 20) -> plt.Figure | None:
    df = df.head(n).iloc[::-1]
    if df.empty:
        return None
    labels = _label_for(df)
    tier_col = df["tier"] if "tier" in df.columns else pd.Series(["low"] * len(df))
    colors = [_tier_color(t) for t in tier_col]
    fig, ax = plt.subplots(figsize=(10, max(6, 0.3 * len(df))))
    bars = ax.barh(labels, df["integrated_score"], color=colors, edgecolor="white")
    for bar, score in zip(bars, df["integrated_score"]):
        ax.text(
            bar.get_width() + 0.01,
            bar.get_y() + bar.get_height() / 2,
            f"{score:.3f}",
            va="center",
            fontsize=9,
        )
    ax.set_xlabel("Integrated Score")
    ax.set_title(f"Top {min(n, len(df))} Neural TF Candidates")
    if not df["integrated_score"].empty:
        ax.set_xlim(0, float(df["integrated_score"].max()) * 1.2)
    return fig


def make_evidence_heatmap(df: pd.DataFrame, n: int = 30) -> plt.Figure | None:
    streams = [s for s in SCORE_STREAMS if s in df.columns]
    if not streams:
        return None
    top = df.head(n)
    heat = top.set_index(_label_for(top))[streams]
    fig, ax = plt.subplots(figsize=(8, max(6, 0.3 * len(heat))))
    sns.heatmap(
        heat.astype(float), annot=True, fmt=".2f", cmap="RdYlGn", center=0.5,
        cbar_kws={"label": "Score"}, ax=ax,
    )
    ax.set_title(f"Top {min(n, len(heat))}: Evidence Stream Scores")
    ax.set_xlabel("Evidence Stream")
    return fig


def make_score_vs_reproducibility(df: pd.DataFrame) -> plt.Figure | None:
    if (
        "reproducibility" not in df.columns
        or "integrated_score" not in df.columns
        or "expression" not in df.columns
    ):
        return None
    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(
        df["reproducibility"], df["integrated_score"],
        c=df["expression"], cmap="viridis", s=80, alpha=0.7, edgecolor="white",
    )
    ax.set_xlabel("Reproducibility (atlases supporting / 3)")
    ax.set_ylabel("Integrated Score")
    ax.set_title("Integrated Score vs Reproducibility (colored by Expression)")
    plt.colorbar(sc, ax=ax, label="Expression Score")
    head = df.head(5)
    nm = _label_for(head)
    for (_, row), name in zip(head.iterrows(), nm):
        ax.annotate(
            str(name),
            (row["reproducibility"], row["integrated_score"]),
            xytext=(5, 5), textcoords="offset points", fontsize=8,
        )
    return fig


def make_expression_vs_specificity(df: pd.DataFrame) -> plt.Figure | None:
    if (
        "expression" not in df.columns
        or "specificity" not in df.columns
        or "tier" not in df.columns
    ):
        return None
    fig, ax = plt.subplots(figsize=(8, 6))
    for tier in ["high", "medium", "low"]:
        sub = df[df["tier"] == tier]
        if sub.empty:
            continue
        ax.scatter(
            sub["expression"], sub["specificity"],
            label=tier, alpha=0.6, s=60, color=_tier_color(tier),
        )
    ax.set_xlabel("Expression Score (max log2FC/5)")
    ax.set_ylabel("Specificity (1 / n_clusters)")
    ax.set_title("Expression vs Specificity by Tier")
    ax.legend()
    return fig


def make_evidence_composition(df: pd.DataFrame, n: int = 15) -> plt.Figure | None:
    streams = [s for s in SCORE_STREAMS
               if s in df.columns and df[s].notna().any()]
    if not streams:
        return None
    top = df.head(n).iloc[::-1]
    labels = _label_for(top)
    fig, ax = plt.subplots(figsize=(10, max(6, 0.4 * len(top))))
    bottom = np.zeros(len(top))
    for s in streams:
        vals = top[s].fillna(0.0).to_numpy(dtype=float)
        ax.barh(labels, vals, left=bottom, label=s, color=STREAM_COLORS[s], edgecolor="white")
        bottom += vals
    ax.set_xlabel("Evidence Score")
    ax.set_title(f"Evidence Composition (Top {min(n, len(top))})")
    ax.legend(loc="lower right")
    return fig


def make_score_by_proof_status(df: pd.DataFrame) -> plt.Figure | None:
    if "proof_status" not in df.columns or "integrated_score" not in df.columns:
        return None
    groups = []
    labels = []
    for s in ["known_rnai_validated", "prior_fstf_not_tested", "novel_candidate"]:
        sub = df[df["proof_status"] == s]["integrated_score"].dropna()
        if not sub.empty:
            groups.append(sub.values)
            labels.append(s.replace("_", " ").title())
    if not groups:
        return None
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.boxplot(groups, patch_artist=True, boxprops=dict(facecolor="#e0e0e0"))
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("Integrated Score")
    ax.set_title("Score Distribution by Proof Status")
    return fig


def make_stream_coverage(df: pd.DataFrame) -> plt.Figure | None:
    if "n_streams" not in df.columns:
        return None
    counts = df["n_streams"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(counts.index, counts.values, color="#4c72b0", edgecolor="white")
    for x, v in zip(counts.index, counts.values):
        ax.text(x, v + 0.5, str(v), ha="center")
    ax.set_xlabel("Number of Supporting Evidence Streams")
    ax.set_ylabel("Number of Candidates")
    ax.set_title("Evidence Stream Coverage")
    return fig


def make_all_vs_neural_scores(run_dir: Path) -> plt.Figure | None:
    all_csv = run_dir / "rank.csv"
    neural_csv = run_dir / "rank_neural.csv"
    if not all_csv.exists() or not neural_csv.exists():
        return None
    a = pd.read_csv(all_csv)
    n = pd.read_csv(neural_csv)
    fig, ax = plt.subplots(figsize=(8, 5))
    if a["integrated_score"].dropna().empty or n["integrated_score"].dropna().empty:
        return None
    bins = np.linspace(0, max(a["integrated_score"].max(),
                              n["integrated_score"].max()), 25)
    ax.hist(a["integrated_score"].dropna(), bins=bins, alpha=0.6,
            label=f"All ({len(a)})", color="#7f7f7f", edgecolor="white")
    ax.hist(n["integrated_score"].dropna(), bins=bins, alpha=0.7,
            label=f"Neural-filtered ({len(n)})", color="#d62728", edgecolor="white")
    ax.set_xlabel("Integrated Score")
    ax.set_ylabel("Candidate Count")
    ax.set_title("Score Distribution: All vs Neural-filtered")
    ax.legend()
    return fig


def make_per_stream_neural_contribution(df: pd.DataFrame) -> plt.Figure | None:
    if "neural_enriched" not in df.columns:
        return None
    streams = [s for s in SCORE_STREAMS
               if s in df.columns and df[s].notna().any()]
    if not streams:
        return None
    enriched = df[df["neural_enriched"].fillna(0.0) > 0]
    other = df[df["neural_enriched"].fillna(0.0) == 0]
    means_e = [enriched[s].mean() if not enriched.empty else 0 for s in streams]
    means_o = [other[s].mean() if not other.empty else 0 for s in streams]
    x = np.arange(len(streams))
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(x - 0.2, means_e, width=0.4, label="Neural-enriched", color="#1f77b4")
    ax.bar(x + 0.2, means_o, width=0.4, label="Not enriched", color="#a0a0a0")
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace("_", "\n") for s in streams], rotation=0, fontsize=8)
    ax.set_ylabel("Mean Score")
    ax.set_title("Per-stream Mean: Neural-enriched vs Not")
    ax.legend()
    return fig


# ---------------------------------------------------------------------------
# PNG-export wrappers (script path). Each calls the corresponding
# `make_*` builder and saves the returned Figure to disk.
# ---------------------------------------------------------------------------

def fig_score_distributions(df: pd.DataFrame, out: Path) -> None:
    f = make_score_distributions(df)
    if f is not None:
        _save(f, out, "1_score_distributions")


def fig_tier_distribution(df: pd.DataFrame, out: Path) -> None:
    f = make_tier_distribution(df)
    if f is not None:
        _save(f, out, "2_tier_distribution")


def fig_proof_status(df: pd.DataFrame, out: Path) -> None:
    f = make_proof_status(df)
    if f is not None:
        _save(f, out, "3_proof_status")


def fig_top_candidates(df: pd.DataFrame, out: Path, n: int = 20) -> None:
    f = make_top_candidates(df, n=n)
    if f is not None:
        _save(f, out, "4_top_candidates")


def fig_evidence_heatmap(df: pd.DataFrame, out: Path, n: int = 30) -> None:
    f = make_evidence_heatmap(df, n=n)
    if f is not None:
        _save(f, out, "5_evidence_heatmap")


def fig_score_vs_reproducibility(df: pd.DataFrame, out: Path) -> None:
    f = make_score_vs_reproducibility(df)
    if f is not None:
        _save(f, out, "6_score_vs_reproducibility")


def fig_expression_vs_specificity(df: pd.DataFrame, out: Path) -> None:
    f = make_expression_vs_specificity(df)
    if f is not None:
        _save(f, out, "7_expression_vs_specificity")


def fig_evidence_composition(df: pd.DataFrame, out: Path, n: int = 15) -> None:
    f = make_evidence_composition(df, n=n)
    if f is not None:
        _save(f, out, "8_evidence_composition")


def fig_score_by_proof_status(df: pd.DataFrame, out: Path) -> None:
    f = make_score_by_proof_status(df)
    if f is not None:
        _save(f, out, "9_score_by_proof_status")


def fig_stream_coverage(df: pd.DataFrame, out: Path) -> None:
    f = make_stream_coverage(df)
    if f is not None:
        _save(f, out, "10_stream_coverage")


def fig_neural_vs_all_score(run_dir: Path, out: Path) -> None:
    f = make_all_vs_neural_scores(run_dir)
    if f is not None:
        _save(f, out, "11_all_vs_neural_scores")


def fig_per_stream_neural_contribution(df: pd.DataFrame, out: Path) -> None:
    f = make_per_stream_neural_contribution(df)
    if f is not None:
        _save(f, out, "12_stream_neural_contribution")


def main() -> None:
    args = parse_args()
    run_dir: Path = args.run
    out_dir: Path = args.out

    rank_csv = run_dir / "rank.csv"
    neural_csv = run_dir / "rank_neural.csv"
    json_path = run_dir / "pipeline_results.json"

    if not rank_csv.exists():
        raise SystemExit(
            f"rank.csv not found in {run_dir}.\n"
            "Run `bioforge neuraltf run` first, or pass --run <dir>."
        )

    print(f"Loading {rank_csv.name}...")
    df = pd.read_csv(rank_csv)
    print(f"  {len(df)} candidates with {len(df.columns)} columns: {list(df.columns)}")

    print(f"Generating figures into {out_dir}/")
    fig_score_distributions(df, out_dir)
    fig_tier_distribution(df, out_dir)
    fig_proof_status(df, out_dir)
    fig_top_candidates(df, out_dir, n=20)
    fig_evidence_heatmap(df, out_dir, n=30)
    fig_score_vs_reproducibility(df, out_dir)
    fig_expression_vs_specificity(df, out_dir)
    fig_evidence_composition(df, out_dir, n=15)
    fig_score_by_proof_status(df, out_dir)
    fig_stream_coverage(df, out_dir)
    fig_neural_vs_all_score(run_dir, out_dir)
    fig_per_stream_neural_contribution(df, out_dir)

    print(f"\nAll figures saved to: {out_dir}")
    print("Generated:")
    for p in sorted(out_dir.glob("*.png")):
        print(f"  {p.name}")


if __name__ == "__main__":
    main()
