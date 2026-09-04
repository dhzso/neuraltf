"""Shared style for NeuralTF publication figures. Single-panel, one graph per image."""
from __future__ import annotations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
RUN = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RES = REPO / "projects" / "NeuralTF" / "results"
FIG = REPO / "projects" / "NeuralTF" / "figures"
SUP = FIG / "supplementary"
FIG.mkdir(parents=True, exist_ok=True)
SUP.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi": 300,
})
plt.rcParams["svg.fonttype"] = "none"

C_A, C_B = "#0072B2", "#E69F00"
C_FIXED, C_CENTERED, C_UNIFORM = "#333333", "#56B4E9", "#CC79A7"
C_NEURAL, C_ALL = "#999999", "#DDDDDD"
C_HL = "#D55E00"
STREAM_COLS = ["expression","specificity","reproducibility","rnai",
               "correlation","neural_enriched","neural_specificity",
               "perez_lineage","perez_influence"]
STREAM_C = {"expression":   "#0072B2",
            "specificity":  "#E69F00",
            "reproducibility": "#009E73",
            "rnai":         "#D55E00",
            "correlation":  "#CC79A7",
            "neural_enriched": "#56B4E9",
            "neural_specificity": "#F0E442",
            "perez_lineage": "#999933",
            "perez_influence": "#6699CC"}   # steel blue — Perez 2025 ANANSE influence
STREAM_L = {"expression":   "Expression",
            "specificity":  "Specificity",
            "reproducibility": "Reproducibility",
            "rnai":         "RNAi",
            "correlation":  "Correlation",
            "neural_enriched": "Neural enriched",
            "neural_specificity": "Neural specificity",
            "perez_lineage": "Perez lineage",
            "perez_influence": "Perez influence"}
# expression=0.2, all 8 others=0.1 (matches EvidenceScorer DEFAULT_WEIGHTS)
W = np.array([0.200, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100])

def _nid(df):
    if "gene_id_v6" in df.columns and "gene_id" not in df.columns:
        df = df.rename(columns={"gene_id_v6": "gene_id"})
    return df


def _csv(path):
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"\n[style.py] Required data file not found:\n  {p}\n"
            "Run the upstream pipeline step that generates this file first.\n"
            "  Full pipeline : python scripts/run.py\n"
            "  Downstream    : python scripts/run_downstream.py"
        )
    df = pd.read_csv(p)
    # Guard: candidate tables must be one row per gene (an upstream
    # annotation-join explosion would silently duplicate every point).
    if "gene_id" in df.columns and not df["gene_id"].is_unique:
        n = df["gene_id"].nunique()
        print(f"[style.py] WARNING: {p.name} has {len(df)} rows but only {n} "
              f"unique gene_id — dropping duplicates (upstream join bug)")
        df = df.drop_duplicates(subset="gene_id", keep="first")
    # Only pad stream columns that the CURRENT pipeline actually defines.
    # (Padding phantom columns renders blank axes and NaN "nan" cells in
    # figures 04/20/31 when a run predates a stream.)
    present_any = any(s in df.columns for s in ["expression", "specificity", "rnai"])
    if present_any:
        base_streams = [s for s in STREAM_COLS if s in df.columns]
        if len(base_streams) < 9:
            print(f"[style.py] NOTE: {p.name} carries only {len(base_streams)} "
                  f"of 9 streams: {base_streams} (run predates a stream?)")
    return df


def load_all():
    return _csv(RUN / "rank.csv")


def load_neural():
    return _csv(RUN / "rank_neural.csv")


def load_top10(f="top10_neural_tfs_prioritized.csv"):
    return _nid(_csv(RES / f))


def load_centered():
    """Load Dirichlet-centered top-10 (5A + 5B)."""
    return _nid(_csv(RES / "dirichlet_centered_top10.csv"))


def load_uniform():
    """Load Dirichlet-uniform top-10 (5A + 5B)."""
    return _nid(_csv(RES / "dirichlet_uniform_top10.csv"))


def load_centered_full():
    """Load Dirichlet-centered full rank (all candidates)."""
    return _nid(_csv(RES / "dirichlet_centered_full_rank.csv"))


def load_uniform_full():
    """Load Dirichlet-uniform full rank (all candidates)."""
    return _nid(_csv(RES / "dirichlet_uniform_full_rank.csv"))


def load_sens_draws():
    return _csv(FIG / "weight_sensitivity_draws.csv")


def load_sens_top10():
    return _csv(FIG / "weight_sensitivity_top10_challengers.csv")


def save(fig, name, dpi=300):
    fig.savefig(FIG / f"{name}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def save_sup(fig, name, dpi=300):
    fig.savefig(SUP / f"{name}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def label(df, gid):
    if "gene_name" in df.columns:
        r = df[df["gene_id"] == gid]
        if len(r) > 0:
            n = r.iloc[0].get("gene_name", "")
            if pd.notna(n) and str(n).strip():
                return str(n)
    return gid
