"""Refined Figure 26: Permutation Null vs Prioritized Top 10 Candidates."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, r"d:/Bioinformatics/projects/NeuralTF/scripts/figures")
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    all_cand = load_all()
    top10 = load_top10()

    stream_cols = [s for s in STREAM_COLS if s in all_cand.columns]
    scores_matrix = all_cand[stream_cols].values.astype(float)
    weights = np.array([W[STREAM_COLS.index(s)] for s in stream_cols])
    weights = weights / weights.sum()

    n_genes = scores_matrix.shape[0]
    n_draws = 50
    rng = np.random.default_rng(42)

    null_scores = []
    for _ in range(n_draws):
        shuffled = scores_matrix.copy()
        for c in range(shuffled.shape[1]):
            col = shuffled[:, c]
            m = ~np.isnan(col)
            idx = np.where(m)[0]
            shuffled[idx, c] = col[idx[rng.permutation(len(idx))]]
        valid = ~np.isnan(shuffled)
        s_fill = np.nan_to_num(shuffled, nan=0.0)
        num = s_fill @ weights
        den = valid.astype(float) @ weights
        den = np.where(den > 0, den, 1.0)
        null_scores.append(num / den)

    null_scores = np.concatenate(null_scores)
    n_null = len(null_scores)

    m = top10.merge(all_cand[["gene_id", "integrated_score"]], on="gene_id", how="left")
    sub_a = m[m["track"] == "A"].sort_values("integrated_score", ascending=True)
    sub_b = m[m["track"] == "B"].sort_values("integrated_score", ascending=True)

    a_scores = sub_a["integrated_score"].values
    b_scores = sub_b["integrated_score"].values

    p99 = float(np.percentile(null_scores, 99))
    med_a = float(np.median(a_scores))
    med_b = float(np.median(b_scores))

    # Exceedance fractions
    frac_a = float((null_scores[:, None] < a_scores[None, :]).mean())
    frac_b = float((null_scores[:, None] < b_scores[None, :]).mean())

    fig, ax = plt.subplots(figsize=(7.5, 4.0), dpi=500)

    # 1. Null distribution histogram & subtle outline
    counts, bins, _ = ax.hist(
        null_scores, bins=55, density=True, color="#CFD8DC", alpha=0.65, edgecolor="#B0BEC5", lw=0.5,
        label=f"Within-stream shuffle null (50 draws × 11,696 = {n_null:,} scores)"
    )
    y_max = float(counts.max())

    # 2. 99th Percentile Threshold
    ax.axvline(x=p99, color="#546E7A", lw=1.2, linestyle="--",
               label=f"99th percentile of null (score = {p99:.2f})")
    ax.text(p99 - 0.015, y_max * 0.95, f"99th %tile\n(s = {p99:.2f})",
            color="#37474F", fontsize=6.2, fontweight="bold", ha="right", va="top")

    # 3. Track B Span & Median Line
    min_b, max_b = float(b_scores.min()), float(b_scores.max())
    ax.axvspan(min_b, max_b, color="#D84315", alpha=0.14)
    ax.axvline(x=med_b, color="#D84315", lw=1.8, linestyle="-",
               label=f"Track B (Unscreened) median ({med_b:.2f})")

    # 4. Track A Span & Median Line
    min_a, max_a = float(a_scores.min()), float(a_scores.max())
    ax.axvspan(min_a, max_a, color="#1565C0", alpha=0.14)
    ax.axvline(x=med_a, color="#1565C0", lw=1.8, linestyle="-",
               label=f"Track A (RNAi-screened) median ({med_a:.2f})")

    # 5. Candidate Markers (Rug / Points)
    # Stagger Track A and Track B vertically so labels and points never collide!
    y_track_a = y_max * 1.10
    y_track_b = y_max * 1.01

    # Plot points
    ax.scatter(a_scores, [y_track_a]*len(a_scores), color="#1565C0", s=36, marker="o", edgecolors="white", linewidth=1.0, zorder=6)
    ax.scatter(b_scores, [y_track_b]*len(b_scores), color="#D84315", s=36, marker="s", edgecolors="white", linewidth=1.0, zorder=6)

    # Distinct Track Headers placed cleanly to the left of their points
    ax.text(0.74, y_track_a, "Track A (RNAi-screened):", fontsize=6.0, fontweight="bold", color="#0D47A1", ha="right", va="center")
    ax.text(0.74, y_track_b, "Track B (Unscreened):", fontsize=6.0, fontweight="bold", color="#BF360C", ha="right", va="center")

    # Median values above points
    ax.text(med_a, y_max * 1.18, f"{med_a:.2f}", fontsize=6.2, fontweight="bold", color="#1565C0", ha="center", va="bottom")
    ax.text(med_b, y_max * 1.18, f"{med_b:.2f}", fontsize=6.2, fontweight="bold", color="#D84315", ha="center", va="bottom")

    # Formatting
    ax.set_xlabel("Integrated Evidence Score", fontsize=7.5, fontweight="bold")
    ax.set_ylabel("Probability Density", fontsize=7.5, fontweight="bold")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0, y_max * 1.25)
    ax.tick_params(axis="both", labelsize=6.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend in upper left, clean text only (NO rectangular box covering data)
    ax.legend(fontsize=5.6, loc="upper left", frameon=False)

    # Clean multi-line title and stats placed strictly ABOVE the graph
    fig.text(0.5, 0.96, "Permutation Null Distribution vs Prioritized Top 10 Candidate Scores",
             fontsize=8.5, fontweight="bold", ha="center")
    fig.text(0.5, 0.915,
             f"Within-stream shuffle null (50 draws × 11,696 candidates = {n_null:,} scores; column marginals preserved)",
             fontsize=6.2, ha="center", color="#555555")
    fig.text(0.5, 0.865,
             f"Track A (RNAi-screened) Median = {med_a:.2f} (P < 0.001)   •   "
             f"Track B (Unscreened) Median = {med_b:.2f} (P < 0.001)   •   "
             f"Null 99th %tile = {p99:.2f}",
             fontsize=6.5, fontweight="bold", ha="center", color="#1A1A1A")

    fig.subplots_adjust(left=0.10, right=0.96, top=0.82, bottom=0.14)
    save(fig, "26_permutation_null")
    print("Successfully built refined 26_permutation_null.png")

if __name__ == "__main__":
    build()
