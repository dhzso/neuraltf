"""Figure 34: ANANSE Regulatory Network Target Specificity.
 
Leverages chromatin accessibility and TF binding profiles across planarian cell fates
to evaluate neuron-specific vs total target capacity for master neural regulators.
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import Counter

def build():
    data_path = RES / "ananse_network_full.csv"
    if not data_path.exists():
        data_path = RES / "ananse_top_regulators.csv"
    if not data_path.exists():
        raise FileNotFoundError(f"{data_path} missing")

    df = pd.read_csv(data_path)
    if "is_ananse_tf" in df.columns:
        tfs = df[df["is_ananse_tf"] == True].copy()
    else:
        tfs = df.copy()

    neural_tfs = tfs[tfs["n_targets_neuron"] > 0].copy()
    other_tfs = tfs[tfs["n_targets_neuron"] == 0].copy()

    fig, ax = plt.subplots(figsize=(W_1COL, 3.2))

    # Reference line: 100% neuron-specific targets (y = x) anchored along specific candidates
    ax.plot([0, 520], [0, 520], color="#888888", ls="--", lw=0.9, zorder=1, label="100% neuron-specific ($y = x$)")

    # Non-neural lineage TFs (at y = 0)
    ax.scatter(
        other_tfs["n_targets_total"],
        other_tfs["n_targets_neuron"],
        color="#CCCCCC",
        edgecolor="#888888",
        linewidth=0.5,
        s=28,
        alpha=0.75,
        label=f"Other lineage TFs ($n = {len(other_tfs)}$)",
        zorder=2,
    )

    # Status colors for neural TFs
    status_colors = {
        "tested": C_A,        # Deep navy
        "not_tested": C_B,             # Muted terracotta
        "known_fstf": "#706E65",  # 2026-09-19: unified Known-FSTF color
    }
    status_labels = {
        "tested": "Tested",
        "not_tested": "Not tested",
        "known_fstf": "Known FSTF",
    }

    # Plot neural TFs by status for clean legend grouping
    for status, group in neural_tfs.groupby("proof_status"):
        color = status_colors.get(status, C_A)
        label_text = status_labels.get(status, status)
        ax.scatter(
            group["n_targets_total"],
            group["n_targets_neuron"],
            color=color,
            edgecolor="#222222",
            linewidth=0.8,
            s=55,
            label=label_text,
            zorder=4,
        )

    # Annotate the neural regulators with computed target counts.
    # 2026-09-19: counts are READ FROM ananse_network_full.csv per row
    # (previously hard-coded "418 / 681" style strings that would go
    # stale if the ANANSE scan is re-run). Offsets are layout only.
    annotations = [
        # (gene_id, label, xytext, ha, va)
        ("dd_Smed_v6_10152_0_1", "dd10152", (16, 4), "left", "center"),
        ("dd_Smed_v6_2442_0_1", "dd2442", (18, -12), "left", "top"),
        ("dd_Smed_v6_8820_0_1", "islet1", (-14, 22), "center", "bottom"),
        ("dd_Smed_v6_11150_0_1", "dd11150", (18, 6), "left", "center"),
        ("dd_Smed_v6_18972_0_1", "dd18972", (18, -16), "left", "top"),
    ]

    for gid, name, offset, ha, va in annotations:
        match = neural_tfs[neural_tfs["v6_id"] == gid]
        if len(match) > 0:
            row = match.iloc[0]
            x_val = row["n_targets_total"]
            y_val = row["n_targets_neuron"]
            txt = f"{name}\n({int(y_val)} / {int(x_val)})"
            ax.annotate(
                txt,
                xy=(x_val, y_val),
                xytext=offset,
                textcoords="offset points",
                fontsize=6.0,
                ha=ha,
                va=va,
                arrowprops=dict(arrowstyle="-", color="#555555", lw=0.6),
                zorder=5,
            )

    ax.set_xlabel("Total predicted target genes ($k_{\\mathrm{total}}$)", fontsize=7.0)
    ax.set_ylabel("Neuron target genes ($k_{\\mathrm{neuron}}$)", fontsize=7.0)
    ax.set_title("ANANSE Regulatory Network Target Specificity", fontsize=8.0, pad=6)
    ax.set_xlim(-30, 950)
    ax.set_ylim(-20, 520)

    # Unframed legend text in upper left so no data area is covered
    ax.legend(
        loc="upper left",
        frameon=False,
        fontsize=6.0,
        handletextpad=0.4,
        borderpad=0.5,
    )

    fig.tight_layout()
    save(fig, "34_ananse_regulatory_network")

if __name__ == "__main__":
    build()
