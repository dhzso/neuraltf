"""Pipeline schematic — visual overview of 5 atlases → 9 streams → scoring → prioritization."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

def build():
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.axis("off")
    ax.set_title("NeuralTF prioritization pipeline schematic",
                 fontweight="bold", fontsize=12, pad=10)

    # Stage 1: Atlases — the REAL five integrated datasets
    atlas_box = FancyBboxPatch((0.2, 1.5), 2.0, 2.0, boxstyle="round,pad=0.1",
                                facecolor="#B3D9FF", edgecolor="#0072B2", linewidth=1.5)
    ax.add_patch(atlas_box)
    ax.text(1.2, 3.3, "Input Atlases", ha="center", fontsize=9, fontweight="bold", color="#0072B2")
    atlases = ["Fincher 2018\n(50.6K cells, v4)", "Plass 2018\n(37.5K cells, v6)",
               "Cui 2023\n(55.0K cells, 8 timepts)", "King 2024\n(TF catalog + RNAi)",
               "Perez 2025\n(lineage + ANANSE)"]
    for i, name in enumerate(atlases):
        ax.text(1.2, 3.0 - i * 0.35, name, ha="center", fontsize=6)

    # Arrow 1
    ax.annotate("", xy=(2.7, 2.5), xytext=(2.2, 2.5),
                arrowprops=dict(arrowstyle="->", color="#333333", lw=1.5))

    # Stage 2: Evidence streams
    stream_box = FancyBboxPatch((2.9, 1.0), 2.2, 3.0, boxstyle="round,pad=0.1",
                                 facecolor="#E6F7E6", edgecolor="#009E73", linewidth=1.5)
    ax.add_patch(stream_box)
    ax.text(4.0, 3.8, "9 Evidence Streams", ha="center", fontsize=9, fontweight="bold", color="#009E73")
    for i, (s, c) in enumerate(STREAM_C.items()):
        y = 3.4 - i * 0.28
        ax.plot(3.1, y, "o", color=c, markersize=5)
        ax.text(3.25, y, STREAM_L[s], fontsize=6, va="center")

    # Arrow 2
    ax.annotate("", xy=(5.6, 2.5), xytext=(5.1, 2.5),
                arrowprops=dict(arrowstyle="->", color="#333333", lw=1.5))

    # Stage 3: Scoring
    score_box = FancyBboxPatch((5.8, 1.8), 1.6, 1.4, boxstyle="round,pad=0.1",
                                facecolor="#FFF2CC", edgecolor="#E69F00", linewidth=1.5)
    ax.add_patch(score_box)
    ax.text(6.6, 2.9, "Evidence Scoring", ha="center", fontsize=9, fontweight="bold", color="#E69F00")
    ax.text(6.6, 2.5, "Weighted sum\n(0.2 + 8×0.1)", ha="center", fontsize=7)
    ax.text(6.6, 2.1, "→ Integrated score", ha="center", fontsize=7)

    # Arrow 3
    ax.annotate("", xy=(7.9, 2.5), xytext=(7.4, 2.5),
                arrowprops=dict(arrowstyle="->", color="#333333", lw=1.5))

    # Stage 4: Prioritization
    prior_box = FancyBboxPatch((8.1, 1.5), 1.6, 2.0, boxstyle="round,pad=0.1",
                                facecolor="#FFE6E6", edgecolor="#D55E00", linewidth=1.5)
    ax.add_patch(prior_box)
    ax.text(8.9, 3.3, "3 Methods", ha="center", fontsize=9, fontweight="bold", color="#D55E00")
    ax.text(8.9, 2.95, "Fixed weights", ha="center", fontsize=7, color=C_FIXED)
    ax.text(8.9, 2.7, "Dirichlet k=40", ha="center", fontsize=7, color=C_CENTERED)
    ax.text(8.9, 2.45, "Dirichlet a=1", ha="center", fontsize=7, color=C_UNIFORM)
    ax.text(8.9, 2.15, "+ same bonuses", ha="center", fontsize=6, color="#555555")
    ax.text(8.9, 1.85, "Track A (RNAi) / B (novel)", ha="center", fontsize=6,
            color="#555555")

    # Bottom annotation
    ax.text(5.0, 0.4, "Pipeline: 5 atlases × 9 evidence streams → weighted scoring → "
            "fixed/centered/uniform methods (shared universe, bonuses, gates)",
            ha="center", fontsize=8, fontstyle="italic", color="#555555")

    save(fig, "22_pipeline_schematic")

if __name__ == "__main__":
    build()
