"""Figure 34: ANANSE Regulatory Network Topology and Target Architecture.

Leverages chromatin accessibility and TF binding profiles across planarian cell fates.
Panel a: Neural vs non-neural target capacity across lineage-defining TFs.
Panel b: Target out-degree of master neural regulators annotated by validation status.
Panel c: Key neural downstream targets shared across top regulators.
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import Counter

def build():
    data_path = RES / "ananse_top_regulators.csv"
    if not data_path.exists():
        raise FileNotFoundError(f"{data_path} missing")

    df = pd.read_csv(data_path)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(W_2COL, 2.9),
                                        gridspec_kw={"width_ratios": [1.1, 1.2, 0.9]})

    # --- Panel a: Neural vs total targets across all TFs ---
    neural_tfs = df[df["n_targets_neuron"] > 0].copy()
    other_tfs = df[df["n_targets_neuron"] == 0].copy()

    ax1.scatter(other_tfs["n_targets_total"], [0]*len(other_tfs),
                color="#CCCCCC", s=25, alpha=0.8, label="Other lineage TFs", zorder=3)
    ax1.scatter(neural_tfs["n_targets_total"], neural_tfs["n_targets_neuron"],
                color=C_A, s=35, edgecolor="#111111", lw=0.6,
                label="Neural regulators", zorder=4)

    # Annotate top neural TFs
    for _, r in neural_tfs.iterrows():
        nm = r["gene_name"] if pd.notna(r["gene_name"]) and r["gene_name"] != "" else r["v6_id"].split("_")[3]
        ax1.annotate(nm, (r["n_targets_total"], r["n_targets_neuron"]),
                     xytext=(4, 2), textcoords="offset points",
                     fontsize=6, fontweight="bold", color=C_A)

    ax1.plot([0, 800], [0, 800], color="#999999", ls=":", lw=0.8, label="100% neural")
    ax1.set_xlabel("Total target genes", fontsize=7.5)
    ax1.set_ylabel("Neuron target genes", fontsize=7.5)
    ax1.set_title("Target specificity", fontsize=8, pad=4)
    ax1.set_xlim(-20, 950)
    ax1.set_ylim(-20, 500)
    ax1.legend(loc="upper left", frameon=False, fontsize=6)
    panel_tag(ax1, "a")

    # --- Panel b: Top neural regulators out-degree ---
    neural_sorted = neural_tfs.sort_values("n_targets_neuron", ascending=True)
    y2 = np.arange(len(neural_sorted))
    labels2 = [r["gene_name"] if pd.notna(r["gene_name"]) and r["gene_name"] != "" else r["v6_id"].split("_")[3]
               for _, r in neural_sorted.iterrows()]
    bar_cols = [C_A if r["proof_status"] == "known_rnai_validated"
                else (C_B if r["proof_status"] == "novel_candidate" else C_NEURAL)
                for _, r in neural_sorted.iterrows()]

    bars = ax2.barh(y2, neural_sorted["n_targets_neuron"], color=bar_cols,
                    height=0.55, edgecolor="none")
    for i, (_, r) in enumerate(neural_sorted.iterrows()):
        ax2.text(r["n_targets_neuron"] + 8, i, f"{int(r['n_targets_neuron'])}",
                 va="center", ha="left", fontsize=6.5, color="#222222")

    ax2.set_yticks(y2)
    ax2.set_yticklabels(labels2, fontsize=7)
    ax2.set_xlabel("Neuron target genes", fontsize=7.5)
    ax2.set_title("Neuron out-degree", fontsize=8, pad=4)
    ax2.set_xlim(0, 500)
    panel_tag(ax2, "b")

    # --- Panel c: Shared downstream neural targets ---
    target_counts = Counter()
    for _, r in neural_tfs.iterrows():
        raw_tgts = str(r["top_5_targets"]).split(";")
        for t in raw_tgts:
            clean_t = t.strip()
            if clean_t and clean_t != "nan":
                target_counts[clean_t] += 1

    common_targets = target_counts.most_common(7)
    if common_targets:
        t_names, t_freqs = zip(*reversed(common_targets))
        t_display = [tn if len(tn) <= 20 else tn[:18] + ".." for tn in t_names]
        y3 = np.arange(len(t_names))
        ax3.barh(y3, t_freqs, color=C_A, height=0.55, edgecolor="none")
        for i, cnt in enumerate(t_freqs):
            ax3.text(cnt + 0.1, i, f"{cnt}/{len(neural_tfs)}",
                     va="center", ha="left", fontsize=6, color="#222222")
        ax3.set_yticks(y3)
        ax3.set_yticklabels(t_display, fontsize=6.5)
        ax3.set_xlabel("Regulating TFs", fontsize=7.5)
        ax3.set_title("Core regulated targets", fontsize=8, pad=4)
        ax3.set_xlim(0, len(neural_tfs) + 1.2)
    panel_tag(ax3, "c")

    fig.tight_layout()
    save(fig, "34_ananse_regulatory_network")

if __name__ == "__main__":
    build()
