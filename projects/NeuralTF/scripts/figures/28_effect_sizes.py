"""Effect sizes — Cliff's delta + Hedges' g bar panel with honesty labels.

2026-09-19: consumes the regenerated effect_sizes.json — the top10 rows
now reference the PUBLISHED dual-track shortlist (not the integrated-score
top-10), and the phenotype_confirmed arm rows (nested honest/honest_strict)
are rendered alongside the neural-vs-rest rows.
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import json

# display label, circular(diagnostic) flag — order = render order
LABELS = {
    "top10_vs_rest":                       ("Published Top-10 vs rest (all streams)", True),
    "top10_vs_rest_honest":                ("Published Top-10 vs rest (label-free)", False),
    "top10_vs_rest_honest_strict":         ("Published Top-10 vs rest (strict label-free)", False),
    "neural_vs_non_neural":                ("Neural vs rest (all streams)", True),
    "neural_vs_non_neural_honest":         ("Neural vs rest (label-free)", False),
    "neural_vs_non_neural_honest_strict":  ("Neural vs rest (strict label-free)", False),
}
# nested arms: json parent key -> (sub-key, label, circular)
NESTED = {
    "phenotype_confirmed": [
        ("honest",        "Phenotype-confirmed vs rest (label-free)", False),
        ("honest_strict", "Phenotype-confirmed vs rest (strict)", False),
    ],
}


def build():
    try:
        data_path = RES / "effect_sizes.json"
        with open(data_path) as f:
            data = json.load(f)

        rows = []  # (stats_dict, label, circular)
        for key, (lbl, circ) in LABELS.items():
            if key in data:
                rows.append((data[key], lbl, circ))
        for parent, subs in NESTED.items():
            for sub, lbl, circ in subs:
                if parent in data and isinstance(data[parent], dict) and sub in data[parent]:
                    rows.append((data[parent][sub], lbl, circ))
        if not rows:
            raise ValueError("effect_sizes.json carries no known comparisons")

        fig, ax = plt.subplots(figsize=(W_15COL, 4.6))

        y = np.arange(len(rows))
        deltas = [r[0]["cliffs_delta"] for r in rows]
        gs = [r[0].get("hedges_g", r[0].get("cohens_d", 0.0)) for r in rows]
        pvals = [r[0].get("p_value", None) for r in rows]
        honest_flags = [not r[2] for r in rows]

        width = 0.32
        for i, (d, h) in enumerate(zip(deltas, honest_flags)):
            ax.barh(y[i] + width/2, d, height=width, color=C_A,
                    alpha=1.0 if h else 0.35, edgecolor="none")
        for i, (g, h) in enumerate(zip(gs, honest_flags)):
            ax.barh(y[i] - width/2, g, height=width, color=C_B,
                    alpha=1.0 if h else 0.35, edgecolor="none")
        from matplotlib.patches import Patch
        handles = [Patch(facecolor=C_A, label="Cliff's δ"),
                   Patch(facecolor=C_B, label="Hedges' g"),
                   Patch(facecolor="#AAAAAA", alpha=0.35,
                         label="Diagnostic (circular — not evidence)")]

        for i, (d, g, p, h) in enumerate(zip(deltas, gs, pvals, honest_flags)):
            suffix = "" if h else "  (diagnostic)"
            ax.text(d + 0.03, y[i] + width/2, f"δ = {d:.2f}{suffix}", va="center", ha="left",
                    fontsize=6.0, color=C_A)
            if p is not None:
                if p < 1e-15:
                    p_math = r"P < 10^{-15}"
                else:
                    base, exp = f"{p:.1e}".split("e")
                    p_math = rf"P = {base} \times 10^{{{int(exp)}}}"
                g_txt = f"g = {g:.2f} (${p_math}$)"
            else:
                g_txt = f"g = {g:.2f}"
            ax.text(g + 0.03, y[i] - width/2, g_txt, va="center", ha="left",
                    fontsize=6.0, color=C_B)

        ax.set_yticks(y)
        ytick_labels = [r[1] + (" †circ." if r[2] else "") for r in rows]
        ax.set_yticklabels(ytick_labels, fontsize=6.8)
        ax.invert_yaxis()
        ax.axvline(x=0, color="#555555", lw=0.6)
        ax.axvline(x=0.5, color="#999999", lw=0.6, linestyle=":", label="Reference (0.5)")
        ax.set_xlabel("Effect size (Cliff's δ and Hedges' g)", fontsize=7.0)
        ax.set_title("Effect Sizes Across Candidate Cohorts (label-free contrasts emphasized)",
                     fontsize=8.0, pad=6)
        ax.legend(handles=handles, fontsize=6.2, loc="lower right", frameon=False)
        ax.set_xlim(-0.05, max(max(deltas), max(gs)) * 1.38)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        save(fig, "28_effect_sizes")

    except FileNotFoundError as e:
        print(f"  [SKIP] {__file__}: {e}")
        return

if __name__ == "__main__":
    build()
