"""Top 10 candidates — graphical information panel with scores and evidence dots."""
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

    # Build merged data
    records = []
    for _, row in top10.iterrows():
        gid = row["gene_id"]
        track = row.get("track","")
        nm = label(neural, gid)
        n_row = neural[neural["gene_id"]==gid]
        c_row = centered[centered["gene_id"]==gid]
        u_row = uniform[uniform["gene_id"]==gid]
        integrated = n_row.iloc[0].get("integrated_score", np.nan) if len(n_row)>0 else np.nan
        comp_fixed = row.get("composite_score", np.nan)
        comp_cent = c_row.iloc[0].get("composite_score", np.nan) if len(c_row)>0 else np.nan
        comp_unif = u_row.iloc[0].get("composite_score", np.nan) if len(u_row)>0 else np.nan
        # Evidence streams
        streams = {}
        for s in STREAM_COLS:
            streams[s] = n_row.iloc[0].get(s, 0) if len(n_row)>0 and pd.notna(n_row.iloc[0].get(s, np.nan)) else 0
        records.append({
            "gene_id": gid, "short": nm, "track": track,
            "integrated": integrated,
            "comp_fixed": comp_fixed, "comp_centered": comp_cent, "comp_uniform": comp_unif,
            "proof": str(row.get("proof_status",""))[:20],
            "domain": str(row.get("interpro_domains", row.get("domains_all",""))).split(";")[0].strip()[:20],
            "ortholog": str(row.get("human_ortholog",""))[:12],
            **streams,
        })
    df = pd.DataFrame(records)

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.axis("off")

    n = len(df)
    bar_h = 0.06
    gap = 0.01
    top_margin = 0.92
    row_h = bar_h + gap

    # Column positions (normalized x)
    col_name = 0.0
    col_integrated = 0.14
    col_comp = 0.22
    col_streams = 0.32
    col_proof = 0.70
    col_domain = 0.80
    col_orth = 0.92

    # Headers
    headers = [("Candidate", col_name), ("Integrated", col_integrated),
               ("Composite", col_comp), ("Evidence streams", col_streams),
               ("Proof status", col_proof), ("Domain", col_domain), ("Ortholog", col_orth)]
    for txt, x in headers:
        ax.text(x, top_margin + 0.02, txt, fontsize=7, fontweight="bold", va="top",
                transform=ax.transAxes, color="#333")
    ax.plot([0, 1], [top_margin, top_margin], color="#CCC", lw=0.5, transform=ax.transAxes)

    for i, r in df.iterrows():
        y = top_margin - i * row_h
        tc = C_A if r["track"] == "A" else C_B
        sym = "o" if r["track"] == "A" else "s"

        # Track symbol + name
        ax.plot(col_name + 0.005, y - bar_h/2, sym, color=tc, markersize=6, transform=ax.transAxes)
        ax.text(col_name + 0.025, y - bar_h/2, f'{r["short"]}', fontsize=8, va="center",
                transform=ax.transAxes, fontweight="bold", color=tc)

        # Integrated score bar
        int_val = r["integrated"]
        if pd.notna(int_val):
            ax.add_patch(FancyBboxPatch(
                (col_integrated, y - bar_h), int_val * 0.08, bar_h,
                boxstyle="round,pad=0.001", facecolor=C_CENTERED, alpha=0.7,
                transform=ax.transAxes))
            ax.text(col_integrated + int_val * 0.08 + 0.003, y - bar_h/2,
                    f"{int_val:.3f}", fontsize=6, va="center", transform=ax.transAxes, fontfamily="monospace")

        # Composite score bar
        comp_val = r["comp_fixed"]
        if pd.notna(comp_val):
            ax.add_patch(FancyBboxPatch(
                (col_comp, y - bar_h), comp_val * 0.08, bar_h,
                boxstyle="round,pad=0.001", facecolor=tc, alpha=0.7,
                transform=ax.transAxes))
            ax.text(col_comp + comp_val * 0.08 + 0.003, y - bar_h/2,
                    f"{comp_val:.3f}", fontsize=6, va="center", transform=ax.transAxes, fontfamily="monospace")

        # Evidence stream dots
        for j, s in enumerate(STREAM_COLS):
            v = r[s]
            alpha = 0.2 + 0.8 * v if v > 0 else 0.15
            c = STREAM_C[s] if v > 0 else "#DDD"
            ax.plot(col_streams + j * 0.055, y - bar_h/2, "o", color=c, markersize=4,
                    alpha=alpha, transform=ax.transAxes)

        # Proof status
        proof_color = "#0072B2" if "validated" in r["proof"].lower() else "#666"
        ax.text(col_proof, y - bar_h/2, r["proof"], fontsize=6, va="center",
                transform=ax.transAxes, color=proof_color)

        # Domain
        ax.text(col_domain, y - bar_h/2, r["domain"], fontsize=6, va="center",
                transform=ax.transAxes, fontstyle="italic", color="#444")

        # Ortholog
        ax.text(col_orth, y - bar_h/2, r["ortholog"], fontsize=6, va="center",
                transform=ax.transAxes, color="#444")

    # Stream legend at bottom
    for j, s in enumerate(STREAM_COLS):
        ax.plot(col_streams + j * 0.055, top_margin - n * row_h - 0.03, "o",
                color=STREAM_C[s], markersize=4, transform=ax.transAxes)
        ax.text(col_streams + j * 0.055 + 0.01, top_margin - n * row_h - 0.03,
                STREAM_L[s][:6], fontsize=5, va="center", transform=ax.transAxes, color="#666")

    # Score legend
    ax.add_patch(FancyBboxPatch((col_integrated, top_margin - n * row_h - 0.06), 0.04, 0.012,
                 boxstyle="round,pad=0.001", facecolor=C_CENTERED, alpha=0.7, transform=ax.transAxes))
    ax.text(col_integrated + 0.045, top_margin - n * row_h - 0.054, "Integrated score", fontsize=5,
            va="center", transform=ax.transAxes, color="#666")
    ax.add_patch(FancyBboxPatch((col_comp, top_margin - n * row_h - 0.06), 0.04, 0.012,
                 boxstyle="round,pad=0.001", facecolor=C_A, alpha=0.7, transform=ax.transAxes))
    ax.text(col_comp + 0.045, top_margin - n * row_h - 0.054, "Composite score (fixed)", fontsize=5,
            va="center", transform=ax.transAxes, color="#666")

    ax.set_title("Prioritized candidates — evidence profile", fontweight="bold", pad=12, fontsize=10, loc="left")
    fig.tight_layout(); save(fig, "05_top10_candidate_atlas")

if __name__=="__main__": build()
