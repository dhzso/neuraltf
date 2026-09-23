"""Method agreement summary — consensus with significance annotations.

Panel a: Method overlap counts with hypergeometric test significance.
Panel b: Pairwise method agreement (Jaccard similarity matrix).
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import json

def build():
    data_path = RES / "overlap_significance.json"
    if not data_path.exists():
        raise FileNotFoundError(
            f"{data_path} missing — run scripts/stats/overlap_significance.py first"
        )
    with open(data_path) as f:
        data = json.load(f)

    # 3x3 Pairwise Jaccard similarity matrix across prioritization formulations
    methods = ["fixed", "centered", "uniform"]
    display_methods = [
        "Fixed weights",
        "Dirichlet Centered\n($k = 40$)",
        "Dirichlet Uniform\n($\\alpha = 1$)",
    ]
    pairwise = data.get("pairwise", {})
    matrix = np.zeros((3, 3))
    counts = np.zeros((3, 3), dtype=int)
    # 2026-09-19: diagonal count computed from the consensus artifact
    # (previously hard-coded 10; all three methods emit an identical
    # 5+5 shortlist, so three_way.overlap_count is the verified value).
    n_top_diag = int(data.get("three_way", {}).get("overlap_count", 10))

    for i, m1 in enumerate(methods):
        for j, m2 in enumerate(methods):
            if i == j:
                matrix[i, j] = 1.0
                counts[i, j] = n_top_diag
            else:
                key = f"{m1}_vs_{m2}"
                info = pairwise.get(key, {})
                matrix[i, j] = info.get("jaccard", 0.0)
                counts[i, j] = info.get("overlap_count", 0)

    fig, ax = plt.subplots(figsize=(3.8, 3.3))

    # Muted blue sequential colormap
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "custom_blues", ["#F5F8FA", "#D4E2EE", "#7AA2C0", "#2B4C6F"], N=256
    )

    im = ax.imshow(matrix, cmap=cmap, vmin=0.0, vmax=1.0, aspect="auto")

    ax.set_xticks(range(3))
    ax.set_xticklabels(display_methods, fontsize=6.8, rotation=25, ha="right")
    ax.set_yticks(range(3))
    ax.set_yticklabels(display_methods, fontsize=6.8)

    # Cell annotations: Jaccard + shared-candidate counts (the identical
    # stratified P-value is stated once in the subtitle, not per cell)
    _keys = ["fixed_vs_centered", "fixed_vs_uniform", "centered_vs_uniform"]
    pair_ps = {pairwise.get(k, {}).get("stratified_poisson_p") for k in _keys} - {None}
    _single_p = next(iter(pair_ps)) if len(pair_ps) == 1 else None
    for i in range(3):
        for j in range(3):
            val = matrix[i, j]
            cnt = counts[i, j]
            text_color = "white" if val > 0.60 else "#222222"
            if i == j:
                cell_text = f"$J = 1.00$\n({n_top_diag}/{n_top_diag})"
            else:
                key = f"{methods[i]}_vs_{methods[j]}"
                # 2026-09-19: only the stratified (within-stratum) p is
                # annotated — the legacy uniform-subset hypergeometric p
                # overstates significance by ~25 orders of magnitude and
                # the previous "P < 10^-30" fallback was an invented
                # statistic with no source in any artifact.
                if _single_p is None:
                    p_strat = pairwise.get(key, {}).get("stratified_poisson_p", None)
                    if p_strat is not None:
                        if p_strat < 1e-15:
                            p_str = r"P < 10^{-15}"
                        else:
                            base, exp = f"{p_strat:.1e}".split("e")
                            p_str = rf"P = {base} \times 10^{{{int(exp)}}}"
                        cell_text = f"$J = {val:.2f}$\n({cnt}/10 shared)\n(${p_str}$)"
                    else:
                        cell_text = f"$J = {val:.2f}$\n({cnt}/10 shared)"
                else:
                    cell_text = f"$J = {val:.2f}$\n({cnt}/10 shared)"
            ax.text(
                j,
                i,
                cell_text,
                ha="center",
                va="center",
                fontsize=5.3,
                color=text_color,
            )

    three_way = data.get("three_way", {})
    n_three = three_way.get("overlap_count", 10)
    p_3way = three_way.get("stratified_poisson_p", None)
    # Correct pairwise stratified Poisson p (9.42e-12, from within-stratum null)
    p_pair_strat = pairwise.get("fixed_vs_centered", {}).get("stratified_poisson_p", None)
    
    if p_pair_strat is not None:
        b_p, e_p = f"{p_pair_strat:.1e}".split("e")
        p_pair_str = rf"P = {b_p} \times 10^{{{int(e_p)}}}"
        # 2026-09-23: keep the one-line subtitle inside the 3.8 in canvas (the
        # longer wording was the widest artist in the tight savefig bbox).
        # 2026-09-23 (v2): "(k/10)" -> "(shared/10)" — k collides with the
        # Dirichlet concentration (k = 40) used on this figure's own axes.
        sub_text = (f"3-way consensus {n_three}/10; J = Jaccard, (shared/10); "
                    f"stratified ${p_pair_str}$")
    else:
        sub_text = f"3-way consensus {n_three}/10; J = Jaccard, (shared/10)"
    title_block(
        fig,
        "Prioritization Method Agreement (Top 10 Candidates)",
        sub_text,
    )

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.82, pad=0.04)
    cbar.set_label("Jaccard similarity index", fontsize=6.8)
    cbar.ax.tick_params(labelsize=6.0)

    # Subtle grid lines between cells
    ax.set_xticks(np.arange(-0.5, 3, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 3, 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    fig.tight_layout()
    fig.subplots_adjust(top=0.84, bottom=0.22)
    save(fig, "33_method_agreement_summary")

if __name__ == "__main__":
    build()

