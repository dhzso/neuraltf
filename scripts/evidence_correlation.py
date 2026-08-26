#!/usr/bin/env python
"""Evidence stream correlation analysis for NeuralTF pipeline.

Computes correlation matrix and PCA of the 7 evidence streams
to assess independence assumptions.

Outputs:
- projects/NeuralTF/results/evidence_correlation_heatmap.png
- projects/NeuralTF/results/evidence_pca.png
- projects/NeuralTF/results/evidence_correlation.csv
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

REPO = Path(__file__).resolve().parents[1]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"

STREAMS = ["expression", "specificity", "reproducibility", "rnai",
           "correlation", "neural_enriched", "neural_specificity"]
STREAM_LABELS = {
    "expression": "Expression",
    "specificity": "Specificity",
    "reproducibility": "Reproducibility",
    "rnai": "RNAi",
    "correlation": "Correlation",
    "neural_enriched": "Neural Enriched",
    "neural_specificity": "Neural Specificity",
}


def main():
    rank_path = RUN_DIR / "rank.csv"
    if not rank_path.exists():
        print(f"Error: {rank_path} not found. Run pipeline first.")
        return 1

    rank = pd.read_csv(rank_path)
    print(f"Loaded {len(rank)} candidates from {rank_path}")

    # Extract evidence stream scores
    S = rank[STREAMS].to_numpy(dtype=float)

    # Replace NaN with 0 for correlation (missing = no evidence)
    S_filled = np.where(np.isnan(S), 0.0, S)

    # Correlation matrix
    corr_matrix = pd.DataFrame(S_filled, columns=STREAMS).corr()
    corr_matrix.index = [STREAM_LABELS[s] for s in STREAMS]
    corr_matrix.columns = [STREAM_LABELS[s] for s in STREAMS]

    # Save correlation matrix
    corr_path = RESULTS_DIR / "evidence_correlation.csv"
    corr_matrix.to_csv(corr_path)
    print(f"Correlation matrix saved to {corr_path}")

    # Plot heatmap
    plt.figure(figsize=(8, 6))
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)
    sns.heatmap(corr_matrix, mask=mask, annot=True, fmt=".2f", cmap="RdBu_r",
                center=0, vmin=-1, vmax=1, square=True,
                cbar_kws={"label": "Pearson r"})
    plt.title("Evidence Stream Correlation Matrix\n(Lower triangle shown)")
    plt.tight_layout()
    heatmap_path = RESULTS_DIR / "evidence_correlation_heatmap.png"
    plt.savefig(heatmap_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Heatmap saved to {heatmap_path}")

    # PCA
    # Only use candidates with at least 1 stream
    has_data = ~np.isnan(S).all(axis=1)
    S_pca = S_filled[has_data]

    pca = PCA()
    pca.fit(S_pca)

    # Variance explained
    var_exp = pca.explained_variance_ratio_
    cum_var = np.cumsum(var_exp)

    # PCA scatter (PC1 vs PC2)
    proj = pca.transform(S_pca)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # PC1 vs PC2 scatter
    axes[0].scatter(proj[:, 0], proj[:, 1], alpha=0.5, s=10)
    axes[0].set_xlabel(f"PC1 ({var_exp[0]:.1%} variance)")
    axes[0].set_ylabel(f"PC2 ({var_exp[1]:.1%} variance)")
    axes[0].set_title("PCA of Evidence Streams (All Candidates)")
    axes[0].grid(True, alpha=0.3)

    # Variance explained
    axes[1].bar(range(1, len(var_exp) + 1), var_exp, alpha=0.7, label="Individual")
    axes[1].plot(range(1, len(cum_var) + 1), cum_var, "r-o", label="Cumulative")
    axes[1].set_xlabel("Principal Component")
    axes[1].set_ylabel("Variance Explained")
    axes[1].set_title("PCA Variance Explained")
    axes[1].set_xticks(range(1, len(var_exp) + 1))
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    pca_path = RESULTS_DIR / "evidence_pca.png"
    plt.savefig(pca_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"PCA plot saved to {pca_path}")

    # Loadings (feature contributions to PCs)
    loadings = pd.DataFrame(
        pca.components_.T,
        index=[STREAM_LABELS[s] for s in STREAMS],
        columns=[f"PC{i+1}" for i in range(len(STREAMS))]
    )
    loadings_path = RESULTS_DIR / "evidence_pca_loadings.csv"
    loadings.to_csv(loadings_path)
    print(f"PCA loadings saved to {loadings_path}")

    # Summary stats
    print("\n=== Evidence Stream Correlation Summary ===")
    print(f"Candidates with ≥1 stream: {has_data.sum()} / {len(rank)}")
    print(f"\nCorrelation matrix:")
    print(corr_matrix.round(3).to_string())

    print(f"\nPCA Variance Explained:")
    for i, (v, c) in enumerate(zip(var_exp, cum_var)):
        print(f"  PC{i+1}: {v:.1%} (cumulative: {c:.1%})")

    print(f"\nPC1 Loadings (strongest contributors):")
    pc1_loadings = loadings["PC1"].abs().sort_values(ascending=False)
    for feat, val in pc1_loadings.items():
        print(f"  {feat}: {loadings.loc[feat, 'PC1']:.3f}")

    # Check for high correlations (>0.7)
    print("\n⚠️  High correlations (|r| > 0.7):")
    for i, s1 in enumerate(STREAMS):
        for j, s2 in enumerate(STREAMS):
            if i < j:
                r = corr_matrix.iloc[i, j]
                if abs(r) > 0.7:
                    print(f"  {STREAM_LABELS[s1]} ↔ {STREAM_LABELS[s2]}: r = {r:.3f}")

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())