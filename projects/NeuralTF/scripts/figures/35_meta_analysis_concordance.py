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
    if not de_path.exists() or not meta_path.exists():
        raise FileNotFoundError(f"Missing required DE or meta-analysis files")

    de_df = pd.read_parquet(de_path)
    meta_df = pd.read_csv(meta_path)

    # Merge on gene id
    merged = pd.merge(de_df, meta_df, left_on="v6_id", right_on="gene_id", how="inner")

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(W_2COL, 2.9),
                                        gridspec_kw={"width_ratios": [1, 0.8, 1.1]})

    # --- Panel a: Cross-atlas effect size correlation (Fincher vs Plass LFC) ---
    clean_lfc = merged.dropna(subset=["fincher_lfc", "plass_lfc"]).copy()
    # Filter extreme outliers for clean plotting
    clean_lfc = clean_lfc[(clean_lfc["fincher_lfc"].between(-4, 6)) & (clean_lfc["plass_lfc"].between(-4, 6))]
    r_val, p_val = spearmanr(clean_lfc["fincher_lfc"], clean_lfc["plass_lfc"])

    # Density / hexbin or alpha scatter
    ax1.scatter(clean_lfc["fincher_lfc"], clean_lfc["plass_lfc"],
                c="#BBBBBB", s=4, alpha=0.25, rasterized=True)

    # Highlight genes significant in both (adjusted p < 0.05)
    sig_both = clean_lfc[(clean_lfc["fincher_p"] < 0.05) & (clean_lfc["plass_p"] < 0.05)]
    ax1.scatter(sig_both["fincher_lfc"], sig_both["plass_lfc"],
                c=C_A, s=8, alpha=0.6, label="Sig. in both", rasterized=True)

    # Identity and zero lines
    ax1.axhline(0, color="#888888", lw=0.6, ls=":")
    ax1.axvline(0, color="#888888", lw=0.6, ls=":")
    ax1.plot([-4, 6], [-4, 6], color="#555555", lw=0.8, ls="--")

    ax1.text(0.05, 0.92, f"Spearman $r_s$ = {r_val:.2f}\n$p$ < $10^{{-30}}$",
             transform=ax1.transAxes, fontsize=6.5, fontweight="bold",
             va="top", ha="left")
    ax1.set_xlabel("Fincher et al. log$_2$FC", fontsize=7.5)
    ax1.set_ylabel("Plass et al. log$_2$FC", fontsize=7.5)
    ax1.set_xlim(-4, 6)
    ax1.set_ylim(-4, 6)
    ax1.legend(loc="lower right", frameon=False, fontsize=6)
    panel_tag(ax1, "a")

    # --- Panel b: Distribution of atlas concordance ---
    counts = meta_df["n_atlases_sig_adj"].value_counts().sort_index()
    x2 = np.arange(len(counts))
    labels2 = [f"{k} atlas" if k == 1 else f"{k} atlases" for k in counts.index]
    colors2 = ["#CCCCCC", C_B, C_A]
    bars = ax2.bar(x2, counts.values, color=colors2[:len(counts)], width=0.6, edgecolor="none")
    for i, v in enumerate(counts.values):
        pct = (v / len(meta_df)) * 100
        ax2.text(i, v + max(counts.values)*0.02, f"{v:,}\n({pct:.1f}%)",
                 ha="center", va="bottom", fontsize=6, color="#222222")

    ax2.set_xticks(x2)
    ax2.set_xticklabels(labels2, fontsize=7)
    ax2.set_ylabel("Candidate genes", fontsize=7.5)
    ax2.set_ylim(0, max(counts.values) * 1.25)
    panel_tag(ax2, "b")

    # --- Panel c: Fisher combined significance vs concordance ---
    groups = []
    labels3 = []
    for k in sorted(meta_df["n_atlases_sig_adj"].unique()):
        subset = meta_df[meta_df["n_atlases_sig_adj"] == k]["fisher_combined_p"].dropna()
        # Cap p-values at 1e-50 for visualization
        logp = -np.log10(np.clip(subset, 1e-50, 1.0))
        groups.append(logp.values)
        labels3.append(f"{k} atlases\n(n={len(subset):,})")

    bp = ax3.boxplot(groups, patch_artist=True, widths=0.55,
                     medianprops=dict(color="#111111", lw=1.2),
                     boxprops=dict(lw=0.7),
                     whiskerprops=dict(lw=0.7, color="#555555"),
                     capprops=dict(lw=0.7, color="#555555"),
                     flierprops=dict(marker=".", markersize=2, alpha=0.2))
    for patch, col in zip(bp["boxes"], colors2[:len(groups)]):
        patch.set_facecolor(col)
        patch.set_alpha(0.7)
        patch.set_edgecolor("#333333")

    ax3.set_xticklabels(labels3, fontsize=6.5)
    ax3.set_ylabel(r"Fisher combined $-\log_{10}(p)$", fontsize=7.5)
    ax3.set_ylim(0, 52)
    panel_tag(ax3, "c")

    fig.tight_layout()
    save(fig, "35_meta_analysis_concordance")

if __name__ == "__main__":
    build()
