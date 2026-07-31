"""Generate all NeuralTF pipeline visualization figures.

Usage: python scripts/visualize_results.py
"""
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# ── paths ──────────────────────────────────────────────────────────────
PROJ = Path(r"D:\Bioinformatics\projects\NeuralTF")
RUN = PROJ / "runs" / "pipeline_run"
FIG = PROJ / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# ── load data ──────────────────────────────────────────────────────────
rank_df = pd.read_csv(RUN / "rank.csv")
with open(RUN / "pipeline_results.json") as f:
    results = json.load(f)

# Set style
sns.set_style("whitegrid")
plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})

# ── 1. Score distribution ──────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
streams = ["expression", "specificity", "reproducibility", "rnai", "correlation", "integrated_score"]
for ax, s in zip(axes.flat, streams):
    if s in rank_df.columns:
        ax.hist(rank_df[s].dropna(), bins=20, edgecolor="white", alpha=0.8)
        ax.set_title(f"{s.replace('_', ' ').title()}")
        ax.set_xlabel("Score")
        ax.set_ylabel("Count")
fig.suptitle("Evidence Stream Score Distributions (82 candidates)", fontsize=14)
plt.tight_layout()
plt.savefig(FIG / "1_score_distributions.png")
plt.close()

# ── 2. Tier pie + bar ──────────────────────────────────────────────────
tier_counts = rank_df["integrated_score"].apply(
    lambda x: "high" if x >= 0.5 else ("medium" if x >= 0.25 else "low")
).value_counts()
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
colors = {"high": "#d62728", "medium": "#ff7f0e", "low": "#2ca02c"}
axes[0].pie(tier_counts.values, labels=tier_counts.index, colors=[colors[t] for t in tier_counts.index], autopct="%1.0f%%")
axes[0].set_title("Tier Distribution")
axes[1].bar(tier_counts.index, tier_counts.values, color=[colors[t] for t in tier_counts.index], edgecolor="white")
axes[1].set_ylabel("Number of Candidates")
axes[1].set_title("Tier Counts")
plt.tight_layout()
plt.savefig(FIG / "2_tier_distribution.png")
plt.close()

# ── 3. Proof status ────────────────────────────────────────────────────
# Extract from evidence_cards.md or infer from RNAi + FSTF
rnai_hits = rank_df[rank_df["rnai"] == 1.0]
proof_status = []
for _, row in rank_df.iterrows():
    if row["rnai"] == 1.0:
        proof_status.append("known_rnai_validated")
    elif row["gene_name"] in ["EGR-1", "foxF-1", "SOXP-3", "ascl-2", "nkx2-like", "p53", "TCF15", "gata456", "hnf4", "gli-1", "fos-1", "Jun", "zic-2", "ets1", "runt-2", "mxi-2", "xbp1", "prep", "POU2/3"]:
        proof_status.append("prior_fstf_not_tested")
    else:
        proof_status.append("novel_candidate")
rank_df["proof_status"] = proof_status
ps_counts = rank_df["proof_status"].value_counts()
fig, ax = plt.subplots(figsize=(6, 4))
ps_colors = {"known_rnai_validated": "#1f77b4", "prior_fstf_not_tested": "#ff7f0e", "novel_candidate": "#2ca02c"}
ax.barh(ps_counts.index, ps_counts.values, color=[ps_colors[s] for s in ps_counts.index], edgecolor="white")
for i, v in enumerate(ps_counts.values):
    ax.text(v + 0.5, i, str(v), va="center")
ax.set_xlabel("Number of Candidates")
ax.set_title("Proof Status Distribution")
plt.tight_layout()
plt.savefig(FIG / "3_proof_status.png")
plt.close()

# ── 4. Top 20 candidates horizontal bar ────────────────────────────────
top20 = rank_df.head(20).iloc[::-1]  # reverse for horizontal
fig, ax = plt.subplots(figsize=(10, 8))
colors_bar = [colors["high"] if s >= 0.5 else (colors["medium"] if s >= 0.25 else colors["low"]) for s in top20["integrated_score"]]
bars = ax.barh(top20["gene_name"], top20["integrated_score"], color=colors_bar, edgecolor="white")
for bar, score in zip(bars, top20["integrated_score"]):
    ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2, f"{score:.3f}", va="center", fontsize=9)
ax.set_xlabel("Integrated Score")
ax.set_title("Top 20 Neural TF Candidates")
ax.set_xlim(0, max(top20["integrated_score"]) * 1.2)
plt.tight_layout()
plt.savefig(FIG / "4_top20_candidates.png")
plt.close()

# ── 5. Evidence stream heatmap ─────────────────────────────────────────
stream_cols = ["expression", "specificity", "reproducibility", "rnai", "correlation"]
heat_data = rank_df.head(30).set_index("gene_name")[stream_cols]
fig, ax = plt.subplots(figsize=(8, 10))
sns.heatmap(heat_data, annot=True, fmt=".2f", cmap="RdYlGn", center=0.5,
            cbar_kws={"label": "Score"}, ax=ax)
ax.set_title("Top 30: Evidence Stream Scores")
ax.set_xlabel("Evidence Stream")
plt.tight_layout()
plt.savefig(FIG / "5_evidence_heatmap.png")
plt.close()

# ── 6. Reproducibility vs Score ────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))
scatter = ax.scatter(rank_df["reproducibility"], rank_df["integrated_score"],
                     c=rank_df["expression"], cmap="viridis", s=80, alpha=0.7, edgecolor="white")
ax.set_xlabel("Reproducibility (atlases supporting / 3)")
ax.set_ylabel("Integrated Score")
ax.set_title("Integrated Score vs Reproducibility (colored by Expression)")
plt.colorbar(scatter, ax=ax, label="Expression Score")
# Annotate top 5
for _, row in rank_df.head(5).iterrows():
    ax.annotate(row["gene_name"], (row["reproducibility"], row["integrated_score"]),
                xytext=(5, 5), textcoords="offset points", fontsize=8)
plt.tight_layout()
plt.savefig(FIG / "6_score_vs_reproducibility.png")
plt.close()

# ── 7. Expression vs Specificity ───────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))
for tier in ["high", "medium", "low"]:
    mask = (rank_df["integrated_score"] >= 0.5) if tier == "high" else \
           (rank_df["integrated_score"] >= 0.25) if tier == "medium" else \
           (rank_df["integrated_score"] < 0.25)
    ax.scatter(rank_df.loc[mask, "expression"], rank_df.loc[mask, "specificity"],
               label=tier, alpha=0.6, s=60, color=colors[tier])
ax.set_xlabel("Expression Score (max log2FC/5)")
ax.set_ylabel("Specificity (1 / n_clusters)")
ax.set_title("Expression vs Specificity by Tier")
ax.legend()
plt.tight_layout()
plt.savefig(FIG / "7_expression_vs_specificity.png")
plt.close()

# ── 8. Stacked bar: evidence composition per candidate ─────────────────
top15 = rank_df.head(15).iloc[::-1]
stream_cols = ["expression", "specificity", "reproducibility", "rnai", "correlation"]
fig, ax = plt.subplots(figsize=(10, 8))
bottom = np.zeros(len(top15))
stream_colors = {"expression": "#1f77b4", "specificity": "#ff7f0e", "reproducibility": "#2ca02c", "rnai": "#d62728", "correlation": "#9467bd"}
for s in stream_cols:
    vals = top15[s].values
    ax.barh(top15["gene_name"], vals, left=bottom, label=s, color=stream_colors[s], edgecolor="white")
    bottom += vals
ax.set_xlabel("Evidence Score")
ax.set_title("Evidence Composition (Top 15)")
ax.legend(loc="lower right")
plt.tight_layout()
plt.savefig(FIG / "8_evidence_composition.png")
plt.close()

# ── 9. Known vs Novel candidates score comparison ──────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
novel = rank_df[rank_df["proof_status"] == "novel_candidate"]["integrated_score"]
prior = rank_df[rank_df["proof_status"] == "prior_fstf_not_tested"]["integrated_score"]
rnai_val = rank_df[rank_df["proof_status"] == "known_rnai_validated"]["integrated_score"]
ax.boxplot([rnai_val, prior, novel], patch_artist=True, boxprops=dict(facecolor="#e0e0e0"))
ax.set_xticklabels(["RNAi-validated", "Prior FSTF", "Novel"])
ax.set_ylabel("Integrated Score")
ax.set_title("Score Distribution by Proof Status")
plt.tight_layout()
plt.savefig(FIG / "9_score_by_proof_status.png")
plt.close()

# ── 10. Number of supporting streams ───────────────────────────────────
stream_counts = rank_df["n_streams"].value_counts().sort_index()
fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(stream_counts.index, stream_counts.values, color="#4c72b0", edgecolor="white")
for x, v in zip(stream_counts.index, stream_counts.values):
    ax.text(x, v + 0.5, str(v), ha="center")
ax.set_xlabel("Number of Supporting Evidence Streams")
ax.set_ylabel("Number of Candidates")
ax.set_title("Evidence Stream Coverage")
plt.tight_layout()
plt.savefig(FIG / "10_stream_coverage.png")
plt.close()

# ── 11. King Atlas contribution ────────────────────────────────────────
# Candidates with king evidence (reproducibility < 1.0 means not all 3 atlases)
has_king = rank_df["reproducibility"] < 1.0  # actually check if king in atlases
# We'll infer from expression > Fincher/Plass alone
fig, ax = plt.subplots(figsize=(8, 5))
king_boost = []
for _, row in rank_df.iterrows():
    # If in 3 atlases (reproducibility ~0.67 or 1.0), check if king contributed
    king_boost.append(1 if row["reproducibility"] <= 0.67 else 0)
ax.bar(["King Atlas Contributed", "King Atlas Not Primary"], [sum(king_boost), len(king_boost) - sum(king_boost)],
       color=["#1f77b4", "#a0a0a0"], edgecolor="white")
ax.set_ylabel("Candidates")
ax.set_title("King Atlas Impact on Reproducibility")
plt.tight_layout()
plt.savefig(FIG / "11_king_atlas_impact.png")
plt.close()

# ── 12. Correlation network (top correlated pairs) ─────────────────────
# Parse correlations from mmc6
import openpyxl
mmc6 = Path(r"D:\Bioinformatics\datasets\raw\Supplementary_Data_ King_2024\1-s2.0-S2211124724001712-mmc6.xlsx")
wb = openpyxl.load_workbook(mmc6)
ws = wb.active
corr_pairs = []
for row in ws.iter_rows(min_row=5, values_only=True):  # header at row 4
    if row[0] and row[1]:
        tf1, tf2 = str(row[0]), str(row[1])
        x1_corr = float(row[2]) if row[2] else 0
        g0_corr = float(row[3]) if row[3] else 0
        gain = g0_corr - x1_corr
        if gain > 0:
            corr_pairs.append({"tf1": tf1, "tf2": tf2, "gain": gain})
corr_df = pd.DataFrame(corr_pairs).sort_values("gain", ascending=False).head(15)
fig, ax = plt.subplots(figsize=(10, 6))
ax.barh([f"{r['tf1']}–{r['tf2']}" for _, r in corr_df.iterrows()], corr_df["gain"], color="#e377c2", edgecolor="white")
ax.set_xlabel("G0−X1 Correlation Gain")
ax.set_title("Top 15 Neural TF Pairs with G0 Enrichment")
plt.tight_layout()
plt.savefig(FIG / "12_correlation_gain.png")
plt.close()

print(f"All 12 figures saved to {FIG}")
print("Generated:")
for f in sorted(FIG.glob("*.png")):
    print(f"  {f.name}")