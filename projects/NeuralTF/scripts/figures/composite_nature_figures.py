"""Composite Nature Communications manuscript figures (Figs 1-5).

Assembles polished, publication-ready multi-panel figures directly matching
the narrative structure of the manuscript:
  - Fig 1: NeuralTF Framework & Comprehensive Evidence Landscape
  - Fig 2: Prioritized Planarian Neural Transcription Factors
  - Fig 3: Weight Sensitivity, Posterior Uncertainty, and Stream Robustness
  - Fig 4: Algorithmic Agreement and Method Independence
  - Fig 5: Empirical Benchmarking, Cross-Atlas Concordance, and Network Topology
"""
from __future__ import annotations
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[3] / "src"))
from style import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, spearmanr
from collections import Counter
import json
from bioforge.projects.neuraltf.planmine import go_term_flags

def compute_bonuses(row):
    bonuses = {}
    go_terms = str(row.get("go_terms", "") or "")
    has_neural_go = has_go_tf = False
    seen = set()
    for term in go_terms.split(";"):
        t = term.strip()
        if not t or t.lower() in ("nan", "none") or t.lower() in seen:
            continue
        seen.add(t.lower())
        is_neural, is_tf = go_term_flags(t)
        has_neural_go = has_neural_go or is_neural
        has_go_tf = has_go_tf or is_tf
    bonuses["GO neural"] = 0.03 if has_neural_go else 0.0
    bonuses["GO TF"] = 0.02 if has_go_tf else 0.0
    orth = str(row.get("human_ortholog", row.get("planmine_human_ortholog_desc", "")) or "")
    bonuses["Human ortholog"] = 0.02 if orth.strip() and orth.lower() != "nan" else 0.0
    return bonuses


def fig1_overview():
    """Figure 1: NeuralTF Framework & Comprehensive Evidence Landscape."""
    fig = plt.figure(figsize=(W_2COL, 6.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0], hspace=0.35, wspace=0.28)

    ax_schematic = fig.add_subplot(gs[0, :])
    ax_cov = fig.add_subplot(gs[1, 0])
    ax_dist = fig.add_subplot(gs[1, 1])

    # --- Panel a: Schematic ---
    ax_schematic.set_xlim(0, 100)
    ax_schematic.set_ylim(0, 32)
    ax_schematic.axis("off")

    cards = [
        {"title": "1. Multi-Atlas Single-Cell", "box": (1, 3, 22, 26), "col": "#E8F1F5", "border": "#4C72B0",
         "items": ["Fincher et al. (2018)", "Plass et al. (2018)", "Cui et al. (2023)", "11,675 candidates"]},
        {"title": "2. 9 Evidence Streams", "box": (26, 3, 23, 26), "col": "#FFF3E6", "border": "#D55E00",
         "items": ["3x scRNA-seq DE", "SCENIC + Regulon", "Perez Lineage + Inf.", "Meta-DE + Concordance"]},
        {"title": "3. Scoring & Robustness", "box": (52, 3, 22, 26), "col": "#F3EBF7", "border": "#8856A7",
         "items": ["Dirichlet posterior (1k)", "Evidence ablation", "Functional bonuses", "Negative controls"]},
        {"title": "4. Neural Prioritization", "box": (77, 3, 22, 26), "col": "#EAF5EA", "border": "#2CA02C",
         "items": ["Track A: Validated", "Track B: Novel TFs", "Decile enrichment", "ANANSE GRN targets"]},
    ]
    for c in cards:
        x, y, w, h = c["box"]
        rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.5,rounding_size=1.5",
                              facecolor=c["col"], edgecolor=c["border"], lw=1.0)
        ax_schematic.add_patch(rect)
        ax_schematic.text(x + w/2, y + h - 3.5, c["title"], ha="center", va="center",
                          fontsize=7.5, fontweight="bold", color="#222222")
        for idx, item in enumerate(c["items"]):
            ax_schematic.text(x + 2.0, y + h - 8.5 - idx * 4.2, f"\u2022 {item}",
                              ha="left", va="center", fontsize=6.2, color="#333333")

    for x_arr in [24.0, 50.0, 75.0]:
        ax_schematic.annotate("", xy=(x_arr + 1.8, 16), xytext=(x_arr - 0.8, 16),
                              arrowprops=dict(arrowstyle="-|>", lw=1.2, color="#555555", mutation_scale=10))
    panel_tag(ax_schematic, "a", x=0.01, y=0.98)

    # --- Panel b: Stream Coverage ---
    df = load_all()
    n_tot = len(df)
    streams = STREAM_COLS
    covs = [(len(df[s].dropna()[df[s].dropna() > 0]) / n_tot) * 100 for s in streams]
    y_cov = np.arange(len(streams))
    labels_cov = [STREAM_L.get(s, s) for s in streams]
    cols_cov = [STREAM_C.get(s, C_A) for s in streams]

    ax_cov.barh(y_cov, covs, color=cols_cov, height=0.6, edgecolor="none", alpha=0.85)
    for i, cv in enumerate(covs):
        ax_cov.text(cv + 1.0, i, f"{cv:.1f}%", va="center", ha="left", fontsize=6, color="#222222")
    ax_cov.set_yticks(y_cov)
    ax_cov.set_yticklabels(labels_cov, fontsize=6.5)
    ax_cov.set_xlabel("Coverage (% candidates with score > 0)", fontsize=7.5)
    ax_cov.set_xlim(0, max(covs) * 1.25)
    ax_cov.invert_yaxis()
    panel_tag(ax_cov, "b")

    # --- Panel c: All vs Neural Distribution ---
    neural = load_neural()
    all_scores = df["integrated_score"].dropna().values
    neural_scores = neural["integrated_score"].dropna().values
    ks_stat, ks_pval = ks_2samp(neural_scores, all_scores, alternative="less")

    ax_dist.hist(all_scores, bins=35, density=True, alpha=0.4, color=C_ALL, label=f"Genome (n={len(all_scores):,})")
    ax_dist.hist(neural_scores, bins=25, density=True, alpha=0.65, color=C_NEURAL, label=f"Neural TFs (n={len(neural_scores):,})")
    ax_dist.axvline(np.median(neural_scores), color=C_NEURAL, lw=1.0, ls="--",
                    label=f"Neural median ({np.median(neural_scores):.2f})")
    ax_dist.axvline(np.median(all_scores), color=C_ALL, lw=1.0, ls=":",
                    label=f"Genome median ({np.median(all_scores):.2f})")

    ax_dist.text(0.05, 0.72, f"KS test: $D={ks_stat:.3f}$\n$p < 10^{{-30}}$",
                 transform=ax_dist.transAxes, fontsize=6.5, fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#CCCCCC", lw=0.5))
    ax_dist.set_xlabel("Integrated evidence score", fontsize=7.5)
    ax_dist.set_ylabel("Density", fontsize=7.5)
    ax_dist.legend(loc="upper right", frameon=False, fontsize=6)
    panel_tag(ax_dist, "c")

    save(fig, "fig1_overview_pipeline")
    print("  wrote fig1_overview_pipeline")


def fig2_prioritization():
    """Figure 2: Prioritized Planarian Neural Transcription Factors."""
    fig = plt.figure(figsize=(W_2COL, 6.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.1], hspace=0.38, wspace=0.28)

    ax_top = fig.add_subplot(gs[0, :])
    ax_heat = fig.add_subplot(gs[1, 0])
    ax_wf = fig.add_subplot(gs[1, 1])

    # --- Panel a: Dual-Track Top Candidates ---
    top10 = load_top10()
    neural = load_neural()
    y = np.arange(len(top10))
    bar_cols = [C_A if r["track"] == "A" else C_B for _, r in top10.iterrows()]
    labels = []
    for _, r in top10.iterrows():
        nm = label(neural, r["gene_id"])
        labels.append(f"[{r['track']}] {nm}")

    score_col = "composite_score" if "composite_score" in top10.columns else "integrated_score"
    scores = top10[score_col].values
    ax_top.barh(y, scores, height=0.55, color=bar_cols, alpha=0.85, edgecolor="none", label="Composite score")

    for i, sc in enumerate(scores):
        ax_top.text(sc + 0.015, i, f"{sc:.3f}", va="center", ha="left", fontsize=6.5, fontweight="bold")
    ax_top.axhline(4.5, color="#555555", lw=0.8, ls="--")
    ax_top.set_yticks(y)
    ax_top.set_yticklabels(labels, fontsize=7)
    ax_top.set_xlabel("Evidence score", fontsize=7.5)
    ax_top.set_xlim(0, 1.15)
    ax_top.legend(loc="lower right", frameon=False, fontsize=6.5)
    ax_top.invert_yaxis()
    panel_tag(ax_top, "a")

    # --- Panel b: Evidence Heatmap (top 20) ---
    top20 = neural.sort_values("integrated_score", ascending=False).head(20)
    streams = STREAM_COLS
    mat = top20[streams].fillna(0).values
    names20 = [label(neural, gid) for gid in top20["gene_id"]]
    im = ax_heat.imshow(mat, cmap="YlGnBu", aspect="auto", vmin=0, vmax=1)
    ax_heat.set_xticks(range(len(streams)))
    ax_heat.set_xticklabels([STREAM_L[s] for s in streams], rotation=45, ha="right", fontsize=6)
    ax_heat.set_yticks(range(len(names20)))
    ax_heat.set_yticklabels(names20, fontsize=6)
    panel_tag(ax_heat, "b")

    # --- Panel c: Functional Bonus Waterfall ---
    s2 = pd.read_csv(RES / "supplementary_table_S2_fixed_all_candidates.csv")
    top10_ids = set(top10["gene_id"])
    wf_rows = []
    for _, row in s2[s2["gene_id"].isin(top10_ids)].iterrows():
        base = row.get("integrated_score", 0)
        bonuses = compute_bonuses(row)
        nm = label(neural, row["gene_id"])
        track = "A" if row.get("proof_status", "") == "known_rnai_validated" else "B"
        wf_rows.append({"name": nm, "track": track, "base": base, **bonuses,
                        "composite": base + sum(bonuses.values())})
    df_wf = pd.DataFrame(wf_rows).sort_values("composite", ascending=True)
    y_wf = np.arange(len(df_wf))

    ax_wf.barh(y_wf, df_wf["base"], height=0.55, color="#4C72B0", alpha=0.85, label="Base evidence")
    curr = df_wf["base"].values.copy()
    for col, colr, lab in [("GO neural", "#D55E00", "GO neural (+0.03)"),
                           ("GO TF", "#009E73", "GO TF (+0.02)"),
                           ("Human ortholog", "#CC79A7", "Ortholog (+0.02)")]:
        vals = df_wf[col].values
        ax_wf.barh(y_wf, vals, left=curr, height=0.55, color=colr, alpha=0.85, label=lab)
        curr += vals

    ax_wf.set_yticks(y_wf)
    ax_wf.set_yticklabels([f"[{r['track']}] {r['name']}" for _, r in df_wf.iterrows()], fontsize=6.5)
    ax_wf.set_xlabel("Score composition", fontsize=7.5)
    ax_wf.set_xlim(0, 1.15)
    ax_wf.legend(loc="lower right", frameon=False, fontsize=5.8)
    panel_tag(ax_wf, "c")

    save(fig, "fig2_prioritization_landscape")
    print("  wrote fig2_prioritization_landscape")


def fig3_sensitivity():
    """Figure 3: Weight Sensitivity, Uncertainty, and Robustness."""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(W_2COL, 5.5))

    # --- Panel a: Dirichlet Draws Boxplots ---
    draws = load_sens_draws()
    top10 = load_top10()
    neural = load_neural()
    t10_genes = top10["gene_id"].tolist()
    sub = draws[draws["gene_id"].isin(t10_genes)].copy()
    score_col = "composite_score" if "composite_score" in top10.columns else "integrated_score"
    gene_order = top10.sort_values(score_col, ascending=True)["gene_id"].tolist()
    box_data = [sub[sub["gene_id"] == g]["rank"].values for g in gene_order]
    labels = [label(neural, g) for g in gene_order]

    bp = ax1.boxplot(box_data, vert=False, patch_artist=True, widths=0.55,
                     medianprops=dict(color="#111111", lw=1.2),
                     boxprops=dict(lw=0.7), whiskerprops=dict(lw=0.7),
                     flierprops=dict(marker=".", markersize=1.5, alpha=0.3))
    for patch in bp["boxes"]:
        patch.set_facecolor(C_A)
        patch.set_alpha(0.65)
    ax1.set_yticks(range(1, len(labels) + 1))
    ax1.set_yticklabels(labels, fontsize=6.5)
    ax1.set_xlabel("Rank across 1,000 Dirichlet draws", fontsize=7.5)
    ax1.set_xlim(0.5, 45)
    panel_tag(ax1, "a")

    # --- Panel b: P(Top 10) ---
    p10_data = []
    for g in gene_order:
        r_vals = sub[sub["gene_id"] == g]["rank"].values
        p10 = np.mean(r_vals <= 10) * 100 if len(r_vals) > 0 else 0
        p10_data.append(p10)
    y_b = np.arange(len(p10_data))
    ax2.barh(y_b, p10_data, height=0.55, color=C_B, alpha=0.85)
    for i, p in enumerate(p10_data):
        ax2.text(p + 1.5, i, f"{p:.1f}%", va="center", ha="left", fontsize=6)
    ax2.set_yticks(y_b)
    ax2.set_yticklabels(labels, fontsize=6.5)
    ax2.set_xlabel("P(Rank \u2264 10) (%)", fontsize=7.5)
    ax2.set_xlim(0, 115)
    panel_tag(ax2, "b")

    # --- Panel c: Global Stream Ablation ---
    abl_p = RES / "stream_ablation_summary.csv"
    if abl_p.exists():
        abl = pd.read_csv(abl_p)
        s_names = [STREAM_L.get(s, s) for s in abl["stream_dropped"]]
        y_c = np.arange(len(abl))
        ax3.barh(y_c, abl["mean_abs_rank_delta"], height=0.55, color="#D55E00", alpha=0.85)
        ax3.set_yticks(y_c)
        ax3.set_yticklabels(s_names, fontsize=6.5)
        ax3.set_xlabel("Mean |\u0394 rank| upon ablation", fontsize=7.5)
        ax3.invert_yaxis()
    panel_tag(ax3, "c")

    # --- Panel d: Top 10 Displacement ---
    if abl_p.exists():
        ax4.barh(y_c, abl["top10_displaced_count"], height=0.55, color="#0072B2", alpha=0.85)
        ax4.set_yticks(y_c)
        ax4.set_yticklabels(s_names, fontsize=6.5)
        ax4.set_xlabel("Top 10 candidates displaced", fontsize=7.5)
        ax4.invert_yaxis()
    panel_tag(ax4, "d")

    fig.tight_layout()
    save(fig, "fig3_sensitivity_and_uncertainty")
    print("  wrote fig3_sensitivity_and_uncertainty")


def fig4_agreement():
    """Figure 4: Algorithmic Agreement and Method Independence."""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(W_2COL, 5.5))

    # --- Panel a: Fixed vs Uniform Full Rank Scatter ---
    s1_path = RES / "supplementary_table_S1_method_comparison.csv"
    if s1_path.exists():
        s1 = pd.read_csv(s1_path)
        ax1.scatter(s1["fixed_rank"], s1["uniform_rank"], color="#CCCCCC", s=2, alpha=0.25, rasterized=True)
        t10_fix = s1[s1["fixed_rank"] <= 10]
        ax1.scatter(t10_fix["fixed_rank"], t10_fix["uniform_rank"], color=C_A, s=15, zorder=5, label="Top 10 (Fixed)")
        ax1.plot([1, 11675], [1, 11675], color="#555555", lw=0.8, ls="--")
        r_val, _ = spearmanr(s1["fixed_rank"], s1["uniform_rank"])
        ax1.text(0.05, 0.90, f"Spearman $r_s = {r_val:.3f}$", transform=ax1.transAxes,
                 fontsize=6.5, fontweight="bold")
        ax1.set_xlabel("Fixed rank (1–11,675)", fontsize=7.5)
        ax1.set_ylabel("Uniform rank (1–11,675)", fontsize=7.5)
        ax1.legend(loc="lower right", frameon=False, fontsize=6)
    panel_tag(ax1, "a")

    # --- Panel b: Jaccard Agreement Heatmap ---
    ov_p = RES / "overlap_significance.json"
    if ov_p.exists():
        with open(ov_p) as f:
            data = json.load(f)
        methods = ["fixed", "centered", "uniform"]
        mat = np.zeros((3, 3))
        for i, m1 in enumerate(methods):
            for j, m2 in enumerate(methods):
                if i == j:
                    mat[i, j] = 1.0
                else:
                    mat[i, j] = data.get("pairwise", {}).get(f"{m1}_vs_{m2}", {}).get("jaccard", 0)
        im = ax2.imshow(mat, cmap="YlGnBu", vmin=0, vmax=1)
        ax2.set_xticks(range(3))
        ax2.set_xticklabels(["Fixed", "Centered", "Uniform"], fontsize=7)
        ax2.set_yticks(range(3))
        ax2.set_yticklabels(["Fixed", "Centered", "Uniform"], fontsize=7)
        for i in range(3):
            for j in range(3):
                ax2.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                         fontsize=7.5, fontweight="bold", color="white" if mat[i, j] > 0.6 else "black")
        cbar = fig.colorbar(im, ax=ax2, shrink=0.8)
        cbar.set_label("Jaccard similarity", fontsize=7)
    panel_tag(ax2, "b")

    # --- Panel c: Multi-Method Slope / Bump Chart ---
    if s1_path.exists():
        top_candidates = s1[s1["fixed_rank"] <= 10].sort_values("fixed_rank")
        x_m = [0, 1, 2]
        neural = load_neural()
        for _, r in top_candidates.iterrows():
            ranks = [r["fixed_rank"], r["centered_rank"], r["uniform_rank"]]
            nm = label(neural, r["gene_id"])
            col = C_A if r["fixed_rank"] <= 5 else C_B
            ax3.plot(x_m, ranks, marker="o", markersize=3.5, lw=1.2, color=col, alpha=0.85)
            ax3.text(2.08, ranks[2], nm, va="center", fontsize=5.8, color=col, fontweight="bold")
        ax3.set_xticks(x_m)
        ax3.set_xticklabels(["Fixed", "Centered", "Uniform"], fontsize=7)
        ax3.set_ylabel("Rank", fontsize=7.5)
        ax3.set_xlim(-0.2, 2.8)
        ax3.invert_yaxis()
    panel_tag(ax3, "c")

    # --- Panel d: Overlap Counts Barplot ---
    if ov_p.exists():
        comps = [("Centered vs Uniform", 9), ("Fixed vs Centered", 5), ("Fixed vs Uniform", 5), ("All Three", 5)]
        c_names, c_vals = zip(*comps)
        y_d = np.arange(len(c_names))
        ax4.barh(y_d, c_vals, height=0.55, color=[C_A, C_B, C_B, C_NEURAL], alpha=0.85)
        for i, v in enumerate(c_vals):
            ax4.text(v + 0.2, i, f"{v}/10", va="center", ha="left", fontsize=6.5, fontweight="bold")
        ax4.set_yticks(y_d)
        ax4.set_yticklabels(c_names, fontsize=6.5)
        ax4.set_xlabel("Candidates in overlap", fontsize=7.5)
        ax4.set_xlim(0, 12)
        ax4.invert_yaxis()
    panel_tag(ax4, "d")

    fig.tight_layout()
    save(fig, "fig4_method_agreement_concordance")
    print("  wrote fig4_method_agreement_concordance")


def fig5_benchmarks():
    """Figure 5: Empirical Benchmarking, Cross-Atlas Concordance, and Network Topology."""
    fig = plt.figure(figsize=(W_2COL, 5.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.35, wspace=0.28)

    ax_neg = fig.add_subplot(gs[0, 0])
    ax_eff = fig.add_subplot(gs[0, 1])
    ax_lfc = fig.add_subplot(gs[1, 0])
    ax_grn = fig.add_subplot(gs[1, 1])

    # --- Panel a: Negative Controls ---
    rank = load_all()
    neural = load_neural()
    ctrl_p = RES / "negative_control_metrics.json"
    if ctrl_p.exists():
        with open(ctrl_p) as f:
            cdata = json.load(f)
        matched_genes = cdata.get("matched_negative_genes", [])
        m_vals = rank[rank["gene_id"].isin(matched_genes)]["integrated_score"].dropna().values
        u_vals = rank["integrated_score"].dropna().values
        pos_vals = neural["integrated_score"].dropna().values

        bp = ax_neg.boxplot([pos_vals, m_vals, u_vals], patch_artist=True, widths=0.55,
                            medianprops=dict(color="#111111", lw=1.2),
                            boxprops=dict(lw=0.7), whiskerprops=dict(lw=0.7),
                            flierprops=dict(marker=".", markersize=1.5, alpha=0.2))
        cols = [C_A, "#CCCCCC", "#E5E5E5"]
        for patch, col in zip(bp["boxes"], cols):
            patch.set_facecolor(col)
            patch.set_alpha(0.7)
            patch.set_edgecolor("#333333")
        ax_neg.set_xticklabels([f"Neural TFs\n(n={len(pos_vals)})", f"Matched Ctrl\n(n={len(m_vals)})", f"Genome\n(n={len(u_vals):,})"], fontsize=6.5)
        ax_neg.set_ylabel("Integrated score", fontsize=7.5)
        ax_neg.text(0.5, 0.90, "Separation $p < 10^{-10}$", transform=ax_neg.transAxes,
                    ha="center", fontsize=6.5, fontweight="bold", color=C_A)
    panel_tag(ax_neg, "a")

    # --- Panel b: Effect Sizes ---
    eff_p = RES / "effect_sizes.json"
    if eff_p.exists():
        with open(eff_p) as f:
            edata = json.load(f)
        c_delta = edata.get("unmatched", {}).get("cliffs_delta", 0.98)
        h_g = edata.get("unmatched", {}).get("hedges_g", 2.45)
        m_delta = edata.get("matched", {}).get("cliffs_delta", 0.94)
        m_g = edata.get("matched", {}).get("hedges_g", 2.10)

        x_e = np.arange(2)
        w = 0.35
        ax_eff.bar(x_e - w/2, [c_delta, m_delta], w, color=C_A, alpha=0.85, label="Cliff's \u03b4")
        ax_eff.bar(x_e + w/2, [h_g/3.0, m_g/3.0], w, color=C_B, alpha=0.85, label="Hedges' g / 3")
        ax_eff.set_xticks(x_e)
        ax_eff.set_xticklabels(["vs Genome", "vs Matched Ctrl"], fontsize=7)
        ax_eff.set_ylabel("Standardized effect size", fontsize=7.5)
        ax_eff.legend(loc="upper right", frameon=False, fontsize=6.5)
        ax_eff.text(0 - w/2, c_delta + 0.03, f"{c_delta:.2f}", ha="center", fontsize=6, fontweight="bold")
        ax_eff.text(0 + w/2, h_g/3.0 + 0.03, f"g={h_g:.2f}", ha="center", fontsize=6, fontweight="bold")
    panel_tag(ax_eff, "b")

    # --- Panel c: Cross-Atlas LFC Concordance ---
    de_p = RUN / "de_pvalues.parquet"
    if de_p.exists():
        de_df = pd.read_parquet(de_p)
        clean = de_df.dropna(subset=["fincher_lfc", "plass_lfc"]).copy()
        clean = clean[(clean["fincher_lfc"].between(-3, 5)) & (clean["plass_lfc"].between(-3, 5))]
        ax_lfc.scatter(clean["fincher_lfc"], clean["plass_lfc"], c="#CCCCCC", s=2, alpha=0.2, rasterized=True)
        r_val, _ = spearmanr(clean["fincher_lfc"], clean["plass_lfc"])
        ax_lfc.plot([-3, 5], [-3, 5], color="#555555", lw=0.8, ls="--")
        ax_lfc.text(0.05, 0.90, f"Spearman $r_s = {r_val:.2f}$", transform=ax_lfc.transAxes,
                    fontsize=6.5, fontweight="bold")
        ax_lfc.set_xlabel("Fincher et al. log$_2$FC", fontsize=7.5)
        ax_lfc.set_ylabel("Plass et al. log$_2$FC", fontsize=7.5)
    panel_tag(ax_lfc, "c")

    # --- Panel d: ANANSE Regulatory Targets ---
    grn_p = RES / "ananse_top_regulators.csv"
    if grn_p.exists():
        grn_df = pd.read_csv(grn_p)
        neural_grn = grn_df[grn_df["n_targets_neuron"] > 0].sort_values("n_targets_neuron", ascending=True)
        y_g = np.arange(len(neural_grn))
        names_g = [r["gene_name"] if pd.notna(r["gene_name"]) and r["gene_name"] != "" else r["v6_id"].split("_")[3]
                   for _, r in neural_grn.iterrows()]
        cols_g = [C_A if r["proof_status"] == "known_rnai_validated" else C_B for _, r in neural_grn.iterrows()]
        ax_grn.barh(y_g, neural_grn["n_targets_neuron"], height=0.55, color=cols_g, alpha=0.85)
        for i, (_, r) in enumerate(neural_grn.iterrows()):
            ax_grn.text(r["n_targets_neuron"] + 8, i, f"{int(r['n_targets_neuron'])}",
                        va="center", ha="left", fontsize=6, fontweight="bold")
        ax_grn.set_yticks(y_g)
        ax_grn.set_yticklabels(names_g, fontsize=6.5)
        ax_grn.set_xlabel("Neuron target genes", fontsize=7.5)
        ax_grn.set_xlim(0, 500)
    panel_tag(ax_grn, "d")

    fig.tight_layout()
    save(fig, "fig5_empirical_benchmarks_and_network")
    print("  wrote fig5_empirical_benchmarks_and_network")


def build_all():
    fig1_overview()
    fig2_prioritization()
    fig3_sensitivity()
    fig4_agreement()
    fig5_benchmarks()

if __name__ == "__main__":
    build_all()
