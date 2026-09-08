"""Shared style for NeuralTF publication figures. Single-panel, one graph per image."""
from __future__ import annotations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
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
    "text.color": "#222222",
    "axes.titlesize": 8.5,
    "axes.labelsize": 8,
    "axes.labelcolor": "#222222",
    "axes.linewidth": 0.6,
    "axes.edgecolor": "#333333",
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "axes.titlepad": 5.0,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "xtick.color": "#222222",
    "ytick.color": "#222222",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.fontsize": 7,
    "legend.labelcolor": "#222222",
    "figure.dpi": 500,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
plt.rcParams["svg.fonttype"] = "none"

# Standard Nature Communications column dimensions (inches)
W_1COL = 3.50   # 89 mm
W_15COL = 5.00  # 127 mm
W_2COL = 7.08   # 180 mm

# Muted, classical scientific publication palette
C_A, C_B = "#2B4C6F", "#B04A3E"           # Deep slate navy (benchmark/Track A), Muted terracotta (candidate/Track B)
C_FIXED, C_CENTERED, C_UNIFORM = "#222222", "#467599", "#5F8D76" # Charcoal, Muted steel blue, Muted sage
C_NEURAL, C_ALL = "#687787", "#C8CED6"    # Slate gray, Soft light gray
C_HL = "#9E3D34"                           # Muted crimson accent / reference
STREAM_COLS = ["expression","specificity","reproducibility","rnai",
               "correlation","neural_enriched","neural_specificity",
               "perez_lineage","perez_influence"]
STREAM_C = {"expression":          "#2B4C6F",  # deep navy
            "specificity":         "#4A7C59",  # sage green
            "reproducibility":     "#5C82A6",  # steel blue
            "rnai":                "#B04A3E",  # terracotta brick
            "correlation":         "#7D5A7D",  # muted plum
            "neural_enriched":     "#3A6B7E",  # deep teal
            "neural_specificity":  "#C08A3E",  # muted warm ochre
            "perez_lineage":       "#7C786E",  # warm slate
            "perez_influence":     "#65799B"}  # slate blue
STREAM_L = {"expression":          "Expression",
            "specificity":         "Specificity",
            "reproducibility":     "Reproducibility",
            "rnai":                "RNAi",
            "correlation":         "Correlation",
            "neural_enriched":     "Neural enriched",
            "neural_specificity":  "Neural specificity",
            "perez_lineage":       "Perez lineage",
            "perez_influence":     "Perez influence"}
STREAM_HEATMAP_ORDER = [
    "expression",
    "specificity",
    "neural_enriched",
    "neural_specificity",
    "reproducibility",
    "rnai",
    "perez_lineage",
    "correlation",
    "perez_influence",
]
# expression=0.2, all 8 others=0.1 (matches EvidenceScorer DEFAULT_WEIGHTS)
W = np.array([0.200, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100])

def _nid(df):
    if "gene_id_v6" in df.columns and "gene_id" not in df.columns:
        df = df.rename(columns={"gene_id_v6": "gene_id"})
    return df


def _csv(path, allow_duplicates=False):
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
    # Repeated measurement matrices (e.g. draws) are exempted.
    if not allow_duplicates and "gene_id" in df.columns and not df["gene_id"].is_unique:
        n = df["gene_id"].nunique()
        print(f"[style.py] WARNING: {p.name} has {len(df)} rows but only {n} "
              f"unique gene_id — dropping duplicates (upstream join bug)")
        df = df.drop_duplicates(subset="gene_id", keep="first")
    # Only pad stream columns that the CURRENT pipeline actually defines.
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
    """Load all 1000 weight sensitivity draws without dropping repeated draws per gene."""
    return _csv(FIG / "weight_sensitivity_draws.csv", allow_duplicates=True)


def load_sens_top10():
    return _csv(FIG / "weight_sensitivity_top10_challengers.csv")


def save(fig, name, dpi=500):
    """Save figure in publication-quality 500 DPI PNG."""
    fig.savefig(FIG / f"{name}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def save_sup(fig, name, dpi=500):
    """Save supplementary figure in publication-quality 500 DPI PNG."""
    fig.savefig(SUP / f"{name}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)



def clean_gene_symbol(name: str | None, gid: str | None = None) -> str:
    """Normalize raw GenBank descriptions and contig strings to clean publication gene symbols.

    e.g. 'Schmidtea mediterranea Zeb-1 mRNA, complete cds' -> 'Zeb-1'
         'Schmidtea mediterranea aristaless-like homeobox transcription factor (arx) mRNA, complete cds' -> 'arx'
         'Schmidtea mediterranea snail-1 mRNA, complete cds' -> 'snail-1'
         'soxB1-2 (not deposited)' -> 'soxB1-2'
         'dd_Smed_v6_19255_0_1' -> 'dd19255'
         'dd_Smed_v6_1854_0_1' -> 'dd1854'
    """
    s = str(name).strip() if (name is not None and pd.notna(name)) else ""
    # Strip '(not deposited)'
    s = re.sub(r"\s*\(\s*not\s+deposited\s*\)", "", s, flags=re.IGNORECASE).strip()

    # GenBank description parsing
    if "schmidtea" in s.lower() or "mrna" in s.lower() or "cds" in s.lower():
        # 1. Symbol in parentheses e.g. '... transcription factor (arx) mRNA ...' or '(ascl-2)'
        m_paren = re.search(
            r"\(([^()]+)\)\s*(?:mRNA|cds|protein|partial|complete)", s, flags=re.IGNORECASE
        )
        if m_paren:
            cand = m_paren.group(1).strip()
            if 1 < len(cand) <= 15:
                return cand
        # 2. 'Schmidtea mediterranea <Symbol> mRNA'
        m_gene = re.search(
            r"Schmidtea\s+mediterranea\s+([A-Za-z0-9/_-]+(?:\s+[0-9/_-]+)?)\s+(?:protein\s+)?(?:mRNA|cds)",
            s,
            flags=re.IGNORECASE,
        )
        if m_gene:
            cand = m_gene.group(1).strip()
            return cand
        # 3. '<Symbol> mRNA'
        m_mrna = re.search(
            r"([A-Za-z0-9/_-]+)\s+(?:protein\s+)?mRNA", s, flags=re.IGNORECASE
        )
        if m_mrna:
            cand = m_mrna.group(1).strip()
            if cand.lower() not in ("schmidtea", "mediterranea", "protein", "hypothetical"):
                return cand

    # If s is a valid short biological gene symbol
    if (
        s
        and len(s) < 25
        and not s.lower().startswith("schmidtea")
        and "mrna" not in s.lower()
        and "hypothetical" not in s.lower()
        and s != "-"
    ):
        m_dd = re.search(r"dd_Smed(?:_v[0-9]+)?_([0-9]+)", s)
        if m_dd:
            return f"dd{m_dd.group(1)}"
        return s

    # Fallback to compact contig ID (e.g. dd_Smed_v6_19255_0_1 -> dd19255)
    target_id = str(gid) if gid else s
    m_dd = re.search(r"dd_Smed(?:_v[0-9]+)?_([0-9]+)", target_id)
    if m_dd:
        return f"dd{m_dd.group(1)}"

    return target_id if target_id else "TF"


def label(df, gid):
    raw_name = ""
    if "gene_name" in df.columns:
        r = df[df["gene_id"] == gid]
        if len(r) > 0:
            n = r.iloc[0].get("gene_name", "")
            if pd.notna(n) and str(n).strip() and str(n).strip() != "-":
                raw_name = str(n).strip()
    return clean_gene_symbol(raw_name, gid)


def panel_tag(ax, letter: str, x: float = -0.12, y: float = 1.05, fontsize: float = 8.5):
    """Add a Nature Communications standard bold lowercase panel tag (a, b, c, ...)."""
    ax.text(x, y, letter.lower(), transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="bottom", ha="right")

