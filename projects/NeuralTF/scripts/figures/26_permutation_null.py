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

        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.hist(null_scores, bins=50, density=True, color=C_NEURAL, alpha=0.6, edgecolor="white",
                label="Null: joint row permutation (evidence vectors)")

        if len(real_scores) > 0:
            # Top-10 BY SCORE among genes TESTABLE by the cluster-label
            # permutation (King-table-saturated genes are excluded: their
            # score cannot change under any cluster permutation).
            testable = df
            if "untestable_by_permutation" in df.columns:
                testable = df[~df["untestable_by_permutation"].astype(bool)]
            top_real = np.sort(testable["real_integrated_score"].dropna().values)[-10:] \
                if "real_integrated_score" in testable.columns else np.sort(real_scores)[-10:]
            for rs in top_real:
                ax.axvline(x=rs, color=C_HL, lw=1.2, linestyle="--", alpha=0.8)
            ax.axvline(x=top_real[-1], color=C_HL, lw=1.5, linestyle="--",
                       label="Top-10 testable candidates")

        p_empirical = df["empirical_p"].min() if "empirical_p" in df.columns else (df["empirical_p_shuffled"].min() if "empirical_p_shuffled" in df.columns else 0.001)
        n_untestable = int(df["untestable_by_permutation"].sum()) if "untestable_by_permutation" in df.columns else 0
        n_testable = len(df) - n_untestable
        n_perm = int(df["n_perm"].iloc[0]) if "n_perm" in df.columns else 30
        p_floor = 1.0 / (n_perm + 1)
        ax.text(0.95, 0.95,
                f"min empirical p = {p_empirical:.4f} "
                f"(= floor 1/{n_perm+1} if {p_empirical:.4g} == {p_floor:.4g})\n"
                f"resolution floor at n={n_perm}: {p_floor:.4f}\n"
                f"testable: {n_testable} | table-saturated: {n_untestable}",
                transform=ax.transAxes, ha="right", va="top", fontsize=8,
                fontweight="bold", color=C_HL,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=C_HL, alpha=0.8))

        ax.set_xlabel("Integrated evidence score")
        ax.set_ylabel("Density")
        ax.set_title("Joint row-permutation null vs top cluster-permutation-testable candidates",
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
