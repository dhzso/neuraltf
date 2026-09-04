"""Permutation null distribution with real scores overlaid."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    try:
        p_path = RES / "permutation_pvalues_full.csv"
        if not p_path.exists():
            p_path = RES / "score_shuffling_pvalues.csv"

        if not p_path.exists():
            print(f"  [SKIP] {__file__}: No permutation file found.")
            return

        df = pd.read_csv(p_path)
        real_scores = df["real_integrated_score"].dropna().values if "real_integrated_score" in df.columns else df.get("real_score", pd.Series()).dropna().values
        # Load or generate empirical null distribution
        all_cand = load_all()
        stream_cols = [s for s in STREAM_COLS if s in all_cand.columns]
        scores_matrix = all_cand[stream_cols].values.astype(float)
        weights = W[:len(stream_cols)]
        weights = weights / weights.sum()

        rng = np.random.default_rng(42)
        null_scores = []
        for _ in range(50):
            shuffled = np.apply_along_axis(rng.permutation, 0, scores_matrix)
            valid = ~np.isnan(shuffled)
            s_fill = np.nan_to_num(shuffled, nan=0.0)
            num = s_fill @ weights
            den = valid.astype(float) @ weights
            den = np.where(den > 0, den, 1.0)
            null_scores.extend(num / den)

        null_scores = np.array(null_scores)

        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.hist(null_scores, bins=50, density=True, color=C_NEURAL, alpha=0.6, edgecolor="white",
                label="Null distribution (shuffled streams)")

        if len(real_scores) > 0:
            # Top-10 BY SCORE (not by p-value ordering — the old head(10)
            # on a p-sorted file marked low-scoring genes)
            top_real = np.sort(real_scores)[-10:]
            for rs in top_real:
                ax.axvline(x=rs, color=C_HL, lw=1.2, linestyle="--", alpha=0.8)
            ax.axvline(x=top_real[-1], color=C_HL, lw=1.5, linestyle="--",
                       label="Top-10 observed candidates")

        p_empirical = df["empirical_p"].min() if "empirical_p" in df.columns else (df["empirical_p_shuffled"].min() if "empirical_p_shuffled" in df.columns else 0.001)
        ax.text(0.95, 0.95, f"Empirical p < {p_empirical:.4f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=9,
                fontweight="bold", color=C_HL,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=C_HL, alpha=0.8))

        ax.set_xlabel("Integrated evidence score")
        ax.set_ylabel("Density")
        ax.set_title("Permutation test: top candidate scores significantly exceed empirical null",
                     fontweight="bold", pad=8)
        ax.legend(fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        save(fig, "26_permutation_null")
    except Exception as e:
        print(f"  [ERROR] {__file__}: {e}")
        return

if __name__ == "__main__":
    build()
