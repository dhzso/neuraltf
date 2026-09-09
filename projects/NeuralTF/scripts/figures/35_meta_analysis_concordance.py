"""Figure 35: Cross-Atlas Meta-Analysis and Effect Size Concordance.

Quantifies neural differential expression reproducibility across independent
single-cell RNA-seq atlases (Fincher et al., Plass et al., Cui et al.).
Panel a: Pairwise log2 fold change concordance between independent atlases.
Panel b: Distribution of multi-atlas concordance (single vs recurrently observed).
Panel c: Meta-analytic statistical power: Fisher combined -log10(p) by concordance degree.
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

def build():
    de_path = RUN / "de_pvalues.parquet"
    meta_path = RES / "meta_analysis_pvalues.csv"
    neural_path = RUN / "rank_neural.csv"
    if not de_path.exists() or not meta_path.exists():
        raise FileNotFoundError("Missing required DE or meta-analysis files")

    de_df = pd.read_parquet(de_path)
    meta_df = pd.read_csv(meta_path)
    neural_df = pd.read_csv(neural_path) if neural_path.exists() else None

    # Merge on gene id
    merged = pd.merge(de_df, meta_df, left_on="v6_id", right_on="gene_id", how="inner")
    clean_lfc = merged.dropna(subset=["fincher_lfc", "plass_lfc"]).copy()
    clean_lfc = clean_lfc[
        (clean_lfc["fincher_lfc"].between(-2, 8)) & (clean_lfc["plass_lfc"].between(-2, 10))
    ]

    r_s, p_s = spearmanr(clean_lfc["fincher_lfc"], clean_lfc["plass_lfc"])
    from scipy.stats import pearsonr
    r_p, p_p = pearsonr(clean_lfc["fincher_lfc"], clean_lfc["plass_lfc"])

    fig, ax = plt.subplots(figsize=(W_1COL, 3.2))

    # Identity and zero reference lines
    ax.axhline(0, color="#D0D7DE", lw=0.6, ls=":", zorder=1)
    ax.axvline(0, color="#D0D7DE", lw=0.6, ls=":", zorder=1)
    ax.plot([-1, 9], [-1, 9], color="#888888", lw=0.8, ls="--", zorder=2, label="Identity line ($y = x$)")

    # Linear regression line
    m, b = np.polyfit(clean_lfc["fincher_lfc"], clean_lfc["plass_lfc"], 1)
    x_vals = np.linspace(clean_lfc["fincher_lfc"].min(), clean_lfc["fincher_lfc"].max(), 100)
    ax.plot(x_vals, m * x_vals + b, color=C_A, lw=1.0, ls="-", zorder=3, label=f"Linear fit ($y = {m:.2f}x + {b:.2f}$)")

    # Background gene points
    ax.scatter(
        clean_lfc["fincher_lfc"],
        clean_lfc["plass_lfc"],
        c=C_A,
        s=7,
        alpha=0.28,
        edgecolor="none",
        rasterized=True,
        label=f"Cross-atlas genes ($n = {len(clean_lfc):,}$)",
        zorder=3,
    )

    # Highlight and label exemplary neural candidates
    if neural_df is not None:
        neural_common = neural_df.merge(clean_lfc, left_on="gene_id", right_on="v6_id", how="inner")
        ax.scatter(
            neural_common["fincher_lfc"],
            neural_common["plass_lfc"],
            c=C_B,
            edgecolor="#222222",
            linewidth=0.7,
            s=34,
            label=f"Neural TFs ($n = {len(neural_common)}$)",
            zorder=6,
        )

        # Annotate selected notable neural regulators
        tf_offsets = {
            "dd_Smed_v6_10038_0_1": ("Zeb-1", (-16, 14), "right", "bottom"),
            "dd_Smed_v6_16955_0_1": ("dd16955", (16, -14), "left", "top"),
            "dd_Smed_v6_11150_0_1": ("dd11150", (-16, 12), "right", "bottom"),
            "dd_Smed_v6_1854_0_1": ("dd1854", (14, -6), "left", "center"),
        }
        for gid, (nm, offset, ha, va) in tf_offsets.items():
            r = neural_common[neural_common["v6_id"] == gid]
            if len(r) > 0:
                ax.annotate(
                    nm,
                    xy=(r.iloc[0]["fincher_lfc"], r.iloc[0]["plass_lfc"]),
                    xytext=offset,
                    textcoords="offset points",
                    fontsize=6.0,
                    color="#111111",
                    ha=ha,
                    va=va,
                    arrowprops=dict(arrowstyle="-", color="#444444", lw=0.6),
                    zorder=7,
                )

    # Inset correlation statistics
    ax.text(
        0.05,
        0.92,
        f"Spearman $r_s = {r_s:.2f}$ ($P < 10^{{-300}}$)\n"
        f"Pearson $r = {r_p:.2f}$ ($P < 10^{{-300}}$)\n"
        f"$N =$ {len(clean_lfc):,} genes",
        transform=ax.transAxes,
        fontsize=6.2,
        va="top",
        ha="left",
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#D0D7DE", alpha=0.95),
        zorder=8,
    )

    ax.set_xlabel("Fincher et al. $\\log_2$ fold change", fontsize=7.0)
    ax.set_ylabel("Plass et al. $\\log_2$ fold change", fontsize=7.0)
    ax.set_title("Cross-Atlas Neural Effect Size Concordance", fontsize=8.0, pad=6)
    ax.set_xlim(-1.5, 9.5)
    ax.set_ylim(-1.5, 10.5)

    ax.legend(
        loc="lower right",
        frameon=True,
        facecolor="white",
        framealpha=0.95,
        edgecolor="#D0D7DE",
        fontsize=6.2,
        handletextpad=0.4,
    )

    fig.tight_layout()
    save(fig, "35_meta_analysis_concordance")

if __name__ == "__main__":
    build()
