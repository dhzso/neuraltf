"""Rank conservation across fixed-weight, centered, and uniform Dirichlet methods.

Panel a: Track A benchmarks (RNAi-validated neural TFs) showing near-perfect rank conservation.
Panel b: Track B novel candidates showing dynamic rank trajectories across Dirichlet priors.
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build():
    f10 = load_top10()
    cf = load_centered_full()
    uf = load_uniform_full()
    neural = load_neural()

    # Calculate within-track ranks across methods
    records = []
    for _, r in f10.iterrows():
        gid = r["gene_id"]
        tr = r["track"]
        cf_tr = cf[cf["proof_status"] == r["proof_status"]] if tr == "A" else cf[cf["proof_status"] != "known_rnai_validated"]
        uf_tr = uf[uf["proof_status"] == r["proof_status"]] if tr == "A" else uf[uf["proof_status"] != "known_rnai_validated"]
        c_sub = cf_tr.reset_index(drop=True)
        u_sub = uf_tr.reset_index(drop=True)
        c_rank = c_sub[c_sub["gene_id"] == gid].index[0] + 1 if gid in c_sub["gene_id"].values else np.nan
        u_rank = u_sub[u_sub["gene_id"] == gid].index[0] + 1 if gid in u_sub["gene_id"].values else np.nan
        records.append({
            "gene_id": gid,
            "name": label(neural, gid),
            "track": tr,
            "fixed": r["rank"],
            "centered": c_rank,
            "uniform": u_rank,
        })

    df = pd.DataFrame(records)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_2COL, 3.2), gridspec_kw={"width_ratios": [1, 1.2]})

    methods = ["Fixed", "Centered\n(k=40)", "Uniform\n(\u03b1=1)"]
    x_pos = [0, 1, 2]

    # --- Panel a: Track A ---
    df_a = df[df["track"] == "A"].sort_values("fixed")
    for _, r in df_a.iterrows():
        ranks = [r["fixed"], r["centered"], r["uniform"]]
        ax1.plot(x_pos, ranks, "o-", color=C_A, lw=1.3, markersize=4, alpha=0.85, zorder=3)
        ax1.text(-0.08, ranks[0], f"{r['name']} (#{int(ranks[0])})",
                 ha="right", va="center", fontsize=6.5, color=C_A, fontweight="bold")
        ax1.text(2.08, ranks[2], f"#{int(ranks[2])}",
                 ha="left", va="center", fontsize=6.5, color=C_A, fontweight="bold")

    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(methods, fontsize=7)
    ax1.set_ylabel("Within-track rank", fontsize=7.5)
    ax1.set_ylim(0.5, 7.5)
    ax1.set_yticks(range(1, 8))
    ax1.invert_yaxis()
    ax1.set_xlim(-1.1, 2.7)
    ax1.set_title("Track A (benchmark)", fontsize=8, pad=4)
    panel_tag(ax1, "a")

    # --- Panel b: Track B ---
    df_b = df[df["track"] == "B"].sort_values("fixed")
    for _, r in df_b.iterrows():
        ranks = [r["fixed"], r["centered"], r["uniform"]]
        ax2.plot(x_pos, ranks, "o-", color=C_B, lw=1.3, markersize=4, alpha=0.85, zorder=3)
        disp_name = clean_gene_symbol(r["name"], r["gene_id"])
        ax2.text(-0.08, ranks[0], f"{disp_name} (#{int(ranks[0])})",
                 ha="right", va="center", fontsize=6.5, color=C_B, fontweight="bold")
        ax2.text(2.08, ranks[2], f"#{int(ranks[2])}",
                 ha="left", va="center", fontsize=6.5, color=C_B, fontweight="bold")

    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(methods, fontsize=7)
    ax2.set_ylabel("Within-track rank", fontsize=7.5)
    ax2.set_ylim(0.5, 24.5)
    ax2.set_yticks([1, 5, 10, 15, 20])
    ax2.invert_yaxis()
    ax2.set_xlim(-1.1, 2.7)
    ax2.set_title("Track B (candidate)", fontsize=8, pad=4)
    panel_tag(ax2, "b")

    fig.tight_layout()
    save(fig, "15_method_bumpchart")

if __name__ == "__main__":
    build()


