"""Shared visual style, color mapping, and utility functions for
NeuralTF publication figures.

All main figures import this module to guarantee a single consistent
visual language across the manuscript.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Sequence

# ---------------------------------------------------------------------------
# Paths (portable — no hardcoded Windows paths)
# ---------------------------------------------------------------------------

def _repo_root() -> Path:
    import os
    return Path(os.environ.get("BIOFORGE_REPO_ROOT", Path.cwd()))

RUN_DIR     = _repo_root() / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = _repo_root() / "projects" / "NeuralTF" / "results"
FIG_DIR     = _repo_root() / "projects" / "NeuralTF" / "figures"
SUPP_DIR    = FIG_DIR / "supplementary"

# ---------------------------------------------------------------------------
# Typography — clean sans-serif, Nature/Cell compatible
# ---------------------------------------------------------------------------
FONT = {
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
}
plt.rcParams.update(FONT)
plt.rcParams["svg.fonttype"] = "none"  # editable text in SVG

# ---------------------------------------------------------------------------
# Semantic color palette — colorblind-safe (Okabe-Ito inspired)
# ---------------------------------------------------------------------------

# Track A / B
C_TRACK_A   = "#0072B2"   # strong blue
C_TRACK_B   = "#E69F00"   # amber

# Method colors
C_FIXED     = "#333333"   # near-black (baseline)
C_CENTERED  = "#56B4E9"   # sky blue
C_UNIFORM   = "#CC79A7"   # mauve/pink

# Populations
C_NEURAL    = "#999999"   # mid-gray for neural-filtered
C_ALL249    = "#DDDDDD"   # light-gray for full universe
C_HIGHLIGHT = "#D55E00"   # vermillion for emphasis

# Evidence stream palette (7 distinct, colorblind-safe)
STREAM_COLORS = {
    "expression":        "#0072B2",
    "specificity":       "#E69F00",
    "reproducibility":   "#009E73",
    "rnai":              "#D55E00",
    "correlation":       "#CC79A7",
    "neural_enriched":   "#56B4E9",
    "neural_specificity": "#F0E442",
}
STREAM_ORDER = ["expression", "specificity", "reproducibility",
                "rnai", "correlation", "neural_enriched", "neural_specificity"]
STREAM_LABELS = {
    "expression": "Expression",
    "specificity": "Specificity",
    "reproducibility": "Reproducibility",
    "rnai": "RNAi",
    "correlation": "Correlation",
    "neural_enriched": "Neural enriched",
    "neural_specificity": "Neural specificity",
}

# Proof status
PROOF_COLORS = {
    "validated":       "#0072B2",
    "prior_fstf":      "#56B4E9",
    "prior_fstf_not_tested": "#999999",
    "novel_candidate": "#E69F00",
    "catalog_fstf_not_in_candidates": "#DDDDDD",
}

# Tier
TIER_COLORS = {
    "HIGH":   "#0072B2",
    "MID":    "#E69F00",
    "LOW":    "#999999",
}

# ---------------------------------------------------------------------------
# Figure dimensions (single column = 89 mm ≈ 3.5 in; double = 180 mm ≈ 7.08 in)
# ---------------------------------------------------------------------------
SINGLE_COL = 3.5
DOUBLE_COL = 7.08
SQUARE     = 3.5

def inch(mm: float) -> float:
    return mm / 25.4

# ---------------------------------------------------------------------------
# Style helpers
# ---------------------------------------------------------------------------

def style_ax(ax, title: str = "", xlabel: str = "", ylabel: str = ""):
    """Apply consistent axis styling."""
    ax.set_title(title, fontweight="bold", pad=6, loc="left")
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=2, width=0.5)
    return ax


def panel_letter(ax, letter: str, x: float = -0.12, y: float = 1.05):
    """Add a bold uppercase panel letter (A, B, C …) to an axes."""
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="top", ha="right")


def savefig(fig, name: str, dpi: int = 300, formats: tuple = ("png", "pdf")):
    """Save figure to FIG_DIR in the requested formats."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        path = FIG_DIR / f"{name}.{fmt}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight",
                    facecolor="white", transparent=False)
    plt.close(fig)


def savefig_supp(fig, name: str, dpi: int = 300, formats: tuple = ("png",)):
    """Save supplementary figure to SUPP_DIR."""
    SUPP_DIR.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        path = SUPP_DIR / f"{name}.{fmt}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight",
                    facecolor="white", transparent=False)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def _normalize_gene_id(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize gene_id_v6 → gene_id for consistent downstream use."""
    if "gene_id_v6" in df.columns and "gene_id" not in df.columns:
        df = df.rename(columns={"gene_id_v6": "gene_id"})
    return df


def load_rank_all() -> pd.DataFrame:
    """Load the full 249-candidate rank.csv."""
    return pd.read_csv(RUN_DIR / "rank.csv")


def load_rank_neural() -> pd.DataFrame:
    """Load the 99-candidate rank_neural.csv."""
    return pd.read_csv(RUN_DIR / "rank_neural.csv")


def load_top10_fixed() -> pd.DataFrame:
    return _normalize_gene_id(pd.read_csv(RESULTS_DIR / "top10_neural_tfs_prioritized.csv"))


def load_top10_centered() -> pd.DataFrame:
    return _normalize_gene_id(pd.read_csv(RESULTS_DIR / "dirichlet_top10_prioritized.csv"))


def load_top10_uniform() -> pd.DataFrame:
    return _normalize_gene_id(pd.read_csv(RESULTS_DIR / "dirichlet_uniform_top10.csv"))


def load_uniform_full_rank() -> pd.DataFrame:
    return _normalize_gene_id(pd.read_csv(RESULTS_DIR / "dirichlet_uniform_full_rank.csv"))


def load_uniform_all249_full_rank() -> pd.DataFrame:
    return _normalize_gene_id(pd.read_csv(RESULTS_DIR / "dirichlet_uniform_all249_full_rank.csv"))


def load_weight_sensitivity_draws() -> pd.DataFrame:
    return pd.read_csv(FIG_DIR / "weight_sensitivity_draws.csv")


def load_weight_sensitivity_top10() -> pd.DataFrame:
    return pd.read_csv(FIG_DIR / "weight_sensitivity_top10_challengers.csv")


# ---------------------------------------------------------------------------
# Candidate ordering
# ---------------------------------------------------------------------------

def order_top10(df: pd.DataFrame, score_col: str = "composite_score") -> list[str]:
    """Return gene_id list ordered: Track A first (descending score),
    then Track B (descending score)."""
    a = df[df["track"] == "A"].sort_values(score_col, ascending=False)
    b = df[df["track"] == "B"].sort_values(score_col, ascending=False)
    return a["gene_id"].tolist() + b["gene_id"].tolist()


def short_label(gene_id: str) -> str:
    """dd31217 → dd31217 (already short). Strip 'dd' prefix for display if long."""
    return gene_id


def gene_label(df: pd.DataFrame, gene_id: str) -> str:
    """Return a human-readable label: 'gene_name (gene_id)' or just gene_id."""
    if "gene_name" in df.columns:
        row = df[df["gene_id"] == gene_id]
        if len(row) > 0:
            name = row.iloc[0].get("gene_name", "")
            if pd.notna(name) and str(name).strip():
                return f"{name}"
    return gene_id
