"""Top 10 candidates — comprehensive information panel (one graph, all info)."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from matplotlib.patches import FancyBboxPatch

def build():
    neural = load_neural()
    top10 = load_top10()
    centered = load_centered()
    uniform = load_uniform()
    order = top10.sort_values("composite_score", ascending=False)

    fig, ax = plt.subplots(figsize=(10, 9))
    ax.axis("off")

    headers = ["Candidate", "Track", "Integrated", "Composite", "Centered", "Uniform",
               "Proof", "Domain", "Ortholog", "RNAi"]
    col_x = [0.0, 0.14, 0.22, 0.30, 0.38, 0.46, 0.54, 0.63, 0.78, 0.88]
    for i, h in enumerate(headers):
        ax.text(col_x[i], 0.97, h, fontsize=7, fontweight="bold", va="top",
                transform=ax.transAxes, color="#333")
    ax.plot([0,1],[0.955,0.955], color="#CCC", lw=0.5, transform=ax.transAxes)

    for idx, (_, row) in enumerate(order.iterrows()):
        gid = row["gene_id"]
        y = 0.92 - idx * 0.085
        track = row.get("track","")
        tc = C_A if track=="A" else C_B
        sym = "o" if track=="A" else "s"

        # Name
        nm = label(neural, gid)
        ax.plot(col_x[0]+0.01, y, sym, color=tc, markersize=5, transform=ax.transAxes)
        ax.text(col_x[0]+0.04, y, nm, fontsize=7, va="center", transform=ax.transAxes, fontweight="bold", color=tc)
        # Track
        ax.text(col_x[1], y, f"{'A (RNAi)' if track=='A' else 'B (novel)'}", fontsize=7, va="center", transform=ax.transAxes, color=tc)
        # Scores
        int_s = row.get("integrated_score", np.nan)
        comp_s = row.get("composite_score", np.nan)
        c_row = centered[centered["gene_id"]==gid]
        u_row = uniform[uniform["gene_id"]==gid]
        cent_s = c_row.iloc[0].get("composite_score",np.nan) if len(c_row)>0 else np.nan
        unif_s = u_row.iloc[0].get("composite_score",np.nan) if len(u_row)>0 else np.nan
        for ci, val in [(2, int_s),(3, comp_s),(4, cent_s),(5, unif_s)]:
            txt = f"{val:.3f}" if pd.notna(val) else "—"
            ax.text(col_x[ci], y, txt, fontsize=7, va="center", transform=ax.transAxes,
                    fontfamily="monospace")
        # Proof
        proof = str(row.get("proof_status",""))[:25]
        ax.text(col_x[6], y, proof, fontsize=6, va="center", transform=ax.transAxes,
                color="#0072B2" if "validated" in proof.lower() else "#666")
        # Domain
        dom = str(row.get("interpro_domains", row.get("domains_all","")))
        if dom and dom != "nan": dom = dom.split(";")[0].strip()[:20]
        else: dom = "—"
        ax.text(col_x[7], y, dom, fontsize=6, va="center", transform=ax.transAxes, fontstyle="italic", color="#444")
        # Ortholog
        orth = str(row.get("human_ortholog",""))
        if not orth or orth == "nan": orth = "—"
        ax.text(col_x[8], y, orth[:12], fontsize=6, va="center", transform=ax.transAxes, color="#444")
        # RNAi
        rnai = str(row.get("rnai_phenotype_notes",""))
        if not rnai or rnai == "nan": rnai = "—"
        else: rnai = rnai[:25]
        ax.text(col_x[9], y, rnai, fontsize=5.5, va="center", transform=ax.transAxes, color="#444")

    ax.set_title("Prioritized candidates — integrated evidence profile", fontweight="bold",
                 pad=12, fontsize=10, loc="left")
    fig.tight_layout(); save(fig, "05_top10_candidate_atlas")

if __name__=="__main__": build()
