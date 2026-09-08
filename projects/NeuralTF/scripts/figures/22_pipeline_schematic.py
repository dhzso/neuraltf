"""Pipeline schematic — visual overview of 5 atlases → 9 streams → scoring → prioritization."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

def build():
    fig, ax = plt.subplots(figsize=(W_2COL, 3.8))

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4.8)
    ax.axis("off")

    # Header title
    ax.text(5.0, 4.55, "NeuralTF: Five-Atlas Evidence Synthesis & Target Prioritization Architecture",
            ha="center", va="center", fontsize=9, fontweight="bold", color="#111111")

    # Box 1: Input Atlases
    b1 = FancyBboxPatch((0.2, 0.4), 2.2, 3.8, boxstyle="round,pad=0.08,rounding_size=0.15",
                         facecolor="#F0F7FD", edgecolor="#0072B2", linewidth=0.8)
    ax.add_patch(b1)
    ax.text(1.3, 3.95, "1. Multi-Atlas Inputs", ha="center", fontsize=8, fontweight="bold", color="#0072B2")
    atlases = [
        ("Fincher 2018", "50,562 cells, scRNA-seq"),
        ("Plass 2018", "37,507 cells, scRNA-seq"),
        ("Cui 2023", "55,014 cells, regeneration"),
        ("King 2024", "TF catalog & RNAi screen"),
        ("Perez 2025", "Lineages & ANANSE GRNs")
    ]
    for i, (name, desc) in enumerate(atlases):
        y_a = 3.45 - i * 0.65
        ax.text(0.4, y_a + 0.12, name, fontsize=7, fontweight="bold", color="#222222")
        ax.text(0.4, y_a - 0.10, desc, fontsize=6, color="#555555")

    # Arrow 1 -> 2
    ax.annotate("", xy=(2.75, 2.3), xytext=(2.45, 2.3),
                arrowprops=dict(arrowstyle="->,head_width=0.3,head_length=0.4", color="#555555", lw=1.2))

    # Box 2: Evidence Streams
    b2 = FancyBboxPatch((2.8, 0.4), 2.35, 3.8, boxstyle="round,pad=0.08,rounding_size=0.15",
                         facecolor="#F2FAF6", edgecolor="#009E73", linewidth=0.8)
    ax.add_patch(b2)
    ax.text(3.97, 3.95, "2. Nine Evidence Streams", ha="center", fontsize=8, fontweight="bold", color="#009E73")
    stream_items = [
        ("Expression", "0.20", STREAM_C["expression"]),
        ("Specificity", "0.10", STREAM_C["specificity"]),
        ("Reproducibility", "0.10", STREAM_C["reproducibility"]),
        ("RNAi Phenotype", "0.10", STREAM_C["rnai"]),
        ("Co-expression Gain", "0.10", STREAM_C["correlation"]),
        ("Neural Enriched", "0.10", STREAM_C["neural_enriched"]),
        ("Neural Specificity", "0.10", STREAM_C["neural_specificity"]),
        ("Perez Lineage Class", "0.10", STREAM_C["perez_lineage"]),
        ("ANANSE Influence", "0.10", STREAM_C["perez_influence"])
    ]
    for i, (name, wt, c) in enumerate(stream_items):
        y_s = 3.55 - i * 0.35
        ax.plot(3.05, y_s, "o", color=c, markersize=4.5)
        ax.text(3.22, y_s, name, fontsize=6.5, va="center", color="#222222")
        ax.text(4.95, y_s, f"w={wt}", fontsize=6, va="center", ha="right", color="#666666")

    # Arrow 2 -> 3
    ax.annotate("", xy=(5.50, 2.3), xytext=(5.20, 2.3),
                arrowprops=dict(arrowstyle="->,head_width=0.3,head_length=0.4", color="#555555", lw=1.2))

    # Box 3: Bayesian Scoring & UQ
    b3 = FancyBboxPatch((5.55, 0.4), 2.15, 3.8, boxstyle="round,pad=0.08,rounding_size=0.15",
                         facecolor="#FFF9F0", edgecolor="#E69F00", linewidth=0.8)
    ax.add_patch(b3)
    ax.text(6.62, 3.95, "3. Scoring & UQ", ha="center", fontsize=8, fontweight="bold", color="#E69F00")
    ax.text(6.62, 3.45, r"$S = \frac{\sum w_i s_i}{\sum w_i}$", fontsize=8, ha="center", color="#222222")
    ax.text(6.62, 3.00, "Three Weight Regimes:", fontsize=6.8, ha="center", fontweight="bold", color="#333333")
    ax.text(6.62, 2.65, "• Fixed weights", fontsize=6.5, ha="center", color=C_FIXED)
    ax.text(6.62, 2.35, "• Centered Dirichlet (k=40)", fontsize=6.5, ha="center", color=C_A)
    ax.text(6.62, 2.05, "• Uniform Dirichlet (α=1)", fontsize=6.5, ha="center", color=C_UNIFORM)
    ax.text(6.62, 1.55, "+ Annotation Bonuses:\nGO Neural (+0.03)\nGO TF (+0.02) | Orth (+0.02)",
            fontsize=6, ha="center", color="#555555")

    # Arrow 3 -> 4
    ax.annotate("", xy=(8.05, 2.3), xytext=(7.75, 2.3),
                arrowprops=dict(arrowstyle="->,head_width=0.3,head_length=0.4", color="#555555", lw=1.2))

    # Box 4: Dual-Track Prioritization
    b4 = FancyBboxPatch((8.1, 0.4), 1.7, 3.8, boxstyle="round,pad=0.08,rounding_size=0.15",
                         facecolor="#FBF4FA", edgecolor="#CC79A7", linewidth=0.8)
    ax.add_patch(b4)
    ax.text(8.95, 3.95, "4. Discovery", ha="center", fontsize=8, fontweight="bold", color="#CC79A7")

    # Track A pill
    tA = FancyBboxPatch((8.25, 2.5), 1.4, 1.05, boxstyle="round,pad=0.04,rounding_size=0.1",
                         facecolor=C_A, edgecolor="none", alpha=0.9)
    ax.add_patch(tA)
    ax.text(8.95, 3.25, "TRACK A", ha="center", fontsize=7, fontweight="bold", color="white")
    ax.text(8.95, 2.90, "RNAi-Validated\nBenchmarks", ha="center", fontsize=6, color="white")
    ax.text(8.95, 2.62, "FoxQ2, six6, unc-4...", ha="center", fontsize=5.5, color="#F0F7FD")

    # Track B pill
    tB = FancyBboxPatch((8.25, 1.05), 1.4, 1.05, boxstyle="round,pad=0.04,rounding_size=0.1",
                         facecolor=C_B, edgecolor="none", alpha=0.9)
    ax.add_patch(tB)
    ax.text(8.95, 1.80, "TRACK B", ha="center", fontsize=7, fontweight="bold", color="white")
    ax.text(8.95, 1.45, "Novel Knockout\nCandidates", ha="center", fontsize=6, color="white")
    ax.text(8.95, 1.17, "ptf-4, dd15328...", ha="center", fontsize=5.5, color="#FFF9F0")

    fig.tight_layout()
    save(fig, "22_pipeline_schematic")

if __name__ == "__main__":
    build()

