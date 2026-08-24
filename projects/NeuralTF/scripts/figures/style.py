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
STREAM_COLS = ["expression","specificity","reproducibility","rnai","correlation","neural_enriched","neural_specificity"]
STREAM_C = {"expression":"#0072B2","specificity":"#E69F00","reproducibility":"#009E73",
            "rnai":"#D55E00","correlation":"#CC79A7","neural_enriched":"#56B4E9","neural_specificity":"#F0E442"}
STREAM_L = {"expression":"Expression","specificity":"Specificity","reproducibility":"Reproducibility",
            "rnai":"RNAi","correlation":"Correlation","neural_enriched":"Neural enriched","neural_specificity":"Neural specificity"}
W = np.array([0.211,0.105,0.158,0.158,0.105,0.158,0.105])

def _nid(df):
    if "gene_id_v6" in df.columns and "gene_id" not in df.columns:
        df = df.rename(columns={"gene_id_v6": "gene_id"})
    return df

def _csv(path):
    return pd.read_csv(path)

def load_all():
    return _csv(RUN/"rank.csv")

def load_neural():
    return _csv(RUN/"rank_neural.csv")

def load_top10(f="top10_neural_tfs_prioritized.csv"):
    return _nid(_csv(RES/f))

def load_centered():
    return _nid(_csv(RES/"dirichlet_top10_prioritized.csv"))

def load_uniform():
    return _nid(_csv(RES/"dirichlet_uniform_top10.csv"))

def load_unif99():
    return _nid(_csv(RES/"dirichlet_uniform_full_rank.csv"))

def load_unif249():
    return _nid(_csv(RES/"dirichlet_uniform_all249_full_rank.csv"))

def load_centered99():
    return _nid(_csv(RES/"dirichlet_centered_full_rank.csv"))

def load_sens_draws():
    return _csv(FIG/"weight_sensitivity_draws.csv")

def load_sens_top10():
    return _csv(FIG/"weight_sensitivity_top10_challengers.csv")

def save(fig, name, dpi=300):
    fig.savefig(FIG/f"{name}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)

def save_sup(fig, name, dpi=300):
    fig.savefig(SUP/f"{name}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)

def label(df, gid):
    if "gene_name" in df.columns:
        r = df[df["gene_id"]==gid]
        if len(r)>0:
            n = r.iloc[0].get("gene_name","")
            if pd.notna(n) and str(n).strip(): return str(n)
    return gid
