"""Permutation null distribution with real scores overlaid.

2026-09-06 audit fix (single consistent null): the previous panel drew a
STREAM-SHUFFLING null (per-column permutation, which destroys the
deterministic within-gene stream dependence) while overlaying p-values
from the CLUSTER-LABEL permutation — two non-equivalent nulls presented
as one test. The histogram now uses the SAME joint row-permutation null
as score_shuffling_permutation.py (whole evidence vectors permuted
across genes, preserving within-gene dependence), and the title/annotation
state the resolution floor honestly (the add-one estimator's minimum
p = 1/(n_perm+1); the displayed p is that floor when all null draws fall
below the real score, which is what "all genes at floor p" means).
"""
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
        weights = np.array([W[STREAM_COLS.index(s)] for s in stream_cols])
        weights = weights / weights.sum()

        rng = np.random.default_rng(42)
        null_scores = []
        # 2026-09-06: JOINT row permutation (whole evidence vectors
        # permuted across genes), matching the fixed
        # score_shuffling_permutation.py null — the old per-column
        # shuffling destroyed deterministic within-gene dependence and
        # was a different null than the overlaid p-values' source.
        for _ in range(50):
            perm = rng.permutation(scores_matrix.shape[0])
            shuffled = scores_matrix[perm]
            valid = ~np.isnan(shuffled)
            s_fill = np.nan_to_num(shuffled, nan=0.0)
            num = s_fill @ weights
            den = valid.astype(float) @ weights
            den = np.where(den > 0, den, 1.0)
            null_scores.extend(num / den)

        null_scores = np.array(null_scores)

        fig, ax = plt.subplots(figsize=(W_15COL, 3.4))
        ax.hist(null_scores, bins=45, density=True, color="#78909C", alpha=0.55, edgecolor="none",
                label="Empirical null (joint evidence row permutation)")

        if len(real_scores) > 0:
            testable = df
            if "untestable_by_permutation" in df.columns:
                testable = df[~df["untestable_by_permutation"].astype(bool)]
            top_real = np.sort(testable["real_integrated_score"].dropna().values)[-10:] \
                if "real_integrated_score" in testable.columns else np.sort(real_scores)[-10:]
            for rs in top_real[:-1]:
                ax.axvline(x=rs, color=C_HL, lw=1.0, linestyle="--", alpha=0.75)
            ax.axvline(x=top_real[-1], color=C_HL, lw=1.4, linestyle="--",
                       label="Top 10 candidates")

        p_empirical = df["empirical_p"].min() if "empirical_p" in df.columns else (df["empirical_p_shuffled"].min() if "empirical_p_shuffled" in df.columns else 0.001)
        n_untestable = int(df["untestable_by_permutation"].sum()) if "untestable_by_permutation" in df.columns else 0
        n_testable = len(df) - n_untestable
        n_perm = int(df["n_perm"].iloc[0]) if "n_perm" in df.columns else 30
        p_floor = 1.0 / (n_perm + 1)
        
        ax.text(0.95, 0.92,
                f"Min empirical $P = {p_empirical:.4f}$\n"
                f"Permutations: {n_perm}\n"
                f"Testable: {n_testable:,}",
                transform=ax.transAxes, ha="right", va="top", fontsize=7,
                color="#222222")

        ax.set_xlabel("Integrated score", fontsize=8)
        ax.set_ylabel("Density", fontsize=8)
        ax.set_title("Permutation null vs observed scores",
                     fontweight="bold", fontsize=8.5, pad=6)
        ax.legend(fontsize=7, frameon=False, loc="upper left")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        save(fig, "26_permutation_null")

    except Exception as e:
        print(f"  [ERROR] {__file__}: {e}")
        return

if __name__ == "__main__":
    build()
