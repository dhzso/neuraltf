"""Cross-Species Ortholog Benchmark.

Independent validation of the NeuralTF integrated prioritization scores against
conserved Homo sapiens transcription factors curated from PlanMine BLAST annotations:
- Panel a: Receiver Operating Characteristic (ROC) curve separating conserved
  neural-fate transcription factors (n=47) from non-neural transcription factors (n=54).
- Panel b: Prioritization score distributions comparing neural-fate vs non-neural
  ortholog classes with median annotations and Mann-Whitney U test significance.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_curve, auc


def build():
    gs_path = RES / "ortholog_gold_standard.csv"
    rank_path = RUN / "rank.csv"
    if not gs_path.exists() or not rank_path.exists():
        raise FileNotFoundError("Missing ortholog gold standard or rank.csv")

    gs = pd.read_csv(gs_path)
    rank = pd.read_csv(rank_path).drop_duplicates(subset="gene_id", keep="first")
    m = rank.merge(gs, on="gene_id")

    pos_df = m[m["label"] == "positive"]
    neg_df = m[m["label"] == "negative"]

    pos_scores = pos_df["integrated_score"].to_numpy()
    neg_scores = neg_df["integrated_score"].to_numpy()

    # Binary labels: 1 for positive, 0 for negative
    y_true = np.concatenate([np.ones(len(pos_scores)), np.zeros(len(neg_scores))])
    y_scores = np.concatenate([pos_scores, neg_scores])

    fpr, tpr, _ = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr, tpr)

    stat, pval = mannwhitneyu(pos_scores, neg_scores, alternative="greater")
    pos_med = np.median(pos_scores)
    neg_med = np.median(neg_scores)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_2COL, 3.4), dpi=500)

    # ------------------ PANEL A: ROC Curve ------------------
    panel_tag(ax1, "a", x=-0.14, y=1.04)
    ax1.plot(fpr, tpr, color=C_A, lw=1.8,
             label=f"Neural-fate orthologs (AUC = {roc_auc:.3f})")
    ax1.plot([0, 1], [0, 1], color="#888888", lw=0.8, linestyle=":",
             label="Random classifier (AUC = 0.500)")

    ax1.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=7.0)
    ax1.set_ylabel("True Positive Rate (Sensitivity)", fontsize=7.0)
    ax1.set_title(f"ROC: Ortholog Discrimination (N = {len(m)})", fontsize=7.8, pad=6)
    ax1.legend(loc="lower right", fontsize=5.8, frameon=False)
    ax1.set_xlim([-0.02, 1.02])
    ax1.set_ylim([-0.02, 1.02])
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # ------------------ PANEL B: Score Boxplot & Jitter ------------------
    panel_tag(ax2, "b", x=-0.14, y=1.04)

    rng = np.random.default_rng(42)
    jitter_neg = rng.normal(0, 0.04, size=len(neg_scores))
    jitter_pos = rng.normal(1, 0.04, size=len(pos_scores))

    # Strip plot (individual points)
    ax2.scatter(jitter_neg, neg_scores, color=C_NEURAL, alpha=0.55, s=16, edgecolors="none", zorder=2)
    ax2.scatter(jitter_pos, pos_scores, color=C_A, alpha=0.55, s=16, edgecolors="none", zorder=2)

    # Boxplot
    bp = ax2.boxplot(
        [neg_scores, pos_scores],
        positions=[0, 1],
        widths=0.38,
        patch_artist=True,
        showfliers=False,
        zorder=3,
        boxprops=dict(facecolor="none", edgecolor="#222222", linewidth=1.0),
        whiskerprops=dict(color="#222222", linewidth=1.0),
        capprops=dict(color="#222222", linewidth=1.0),
        medianprops=dict(color=C_B, linewidth=1.8),
    )

    # Median text annotations
    ax2.text(0 + 0.24, neg_med, f"Med={neg_med:.3f}", va="center", ha="left",
             fontsize=6.0, fontweight="bold", color="#444444")
    ax2.text(1 + 0.24, pos_med, f"Med={pos_med:.3f}", va="center", ha="left",
             fontsize=6.0, fontweight="bold", color=C_A)

    # Significance bracket
    y_bar = max(m["integrated_score"].max() + 0.04, 0.98)
    h_tick = 0.02
    ax2.plot([0, 0, 1, 1], [y_bar - h_tick, y_bar, y_bar, y_bar - h_tick], color="#333333", lw=0.9)
    p_str = f"P = {pval:.2e}" if pval < 0.001 else f"P = {pval:.3f}"
    ax2.text(0.5, y_bar + 0.015, f"Mann-Whitney U = {stat:,.0f} ({p_str})",
             ha="center", va="bottom", fontsize=6.0, fontweight="bold", color="#222222")

    ax2.set_xticks([0, 1])
    ax2.set_xticklabels([f"Non-neural\n(n = {len(neg_scores)})", f"Neural-fate\n(n = {len(pos_scores)})"],
                        fontsize=6.8, fontweight="bold")
    ax2.set_ylabel("Integrated Neural TF Score", fontsize=7.0)
    ax2.set_title("Prioritization Score Distribution", fontsize=7.8, pad=6)
    ax2.set_xlim([-0.45, 1.55])
    ax2.set_ylim([m["integrated_score"].min() - 0.05, y_bar + 0.08])
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.suptitle("Cross-Species Ortholog Benchmark: Conserved Human Neural TF Recovery",
                 fontsize=8.5, fontweight="bold", y=0.99)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.14, wspace=0.25)
    save(fig, "37_ortholog_benchmark")


if __name__ == "__main__":
    build()
