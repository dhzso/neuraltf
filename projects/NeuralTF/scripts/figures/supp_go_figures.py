"""GO term figures — comprehensive annotation landscape, namespace, and track comparison.

Produces:
  - fig_s1_go_landscape: Full GO annotation matrix across candidates (Track A vs B).
  - fig_s2_go_namespace_and_track: GO namespace composition and per-track annotation density.
"""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

GO_REF = FIG / "go_term_reference.csv"
GO_MATRIX = SUP / "go_gene_term_matrix_reduced.csv"

def load_go():
    return pd.read_csv(GO_REF), pd.read_csv(GO_MATRIX)

def fig_s1_go_landscape():
    ref, mat = load_go()
    neural = load_neural()
    go_ids = [c for c in mat.columns if c.startswith("GO:")]
    term_neural = dict(zip(ref["go_id"], ref["neural_go"] == "yes"))
    term_tf = dict(zip(ref["go_id"], ref["tf_go"] == "yes"))
    term_name = dict(zip(ref["go_id"], ref["term"]))
    term_count = dict(zip(ref["go_id"], ref["n_of_97_neural_candidates"]))
    go_ids_sorted = sorted(go_ids, key=lambda t: (
        not term_neural.get(t, False), not term_tf.get(t, False),
        -term_count.get(t, 0)))
    go_ids_filtered = [t for t in go_ids_sorted if term_count.get(t, 0) >= 2]

    neural_end = sum(1 for t in go_ids_filtered if term_neural.get(t, False))
    tf_end = neural_end + sum(1 for t in go_ids_filtered[neural_end:] if term_tf.get(t, False))

    proof_map = dict(zip(neural["gene_id"], neural["proof_status"]))
    score_map = dict(zip(neural["gene_id"], neural["integrated_score"]))
    gene_ids = sorted(mat["gene_id"].tolist(),
                      key=lambda g: (proof_map.get(g, "") != "known_rnai_validated",
                                     -score_map.get(g, 0)))
    mat_indexed = mat.set_index("gene_id").loc[gene_ids]
    data = mat_indexed[go_ids_filtered].fillna(0).values

    fig, ax = plt.subplots(figsize=(W_2COL, 9.5))

    cmap = plt.cm.colors.ListedColormap(["#F7F7F7", C_A])
    ax.imshow(data, cmap=cmap, aspect="auto", interpolation="nearest", vmin=0, vmax=1)

    n_genes = len(gene_ids)
    if 0 < neural_end < len(go_ids_filtered):
        ax.plot([neural_end - 0.5, neural_end - 0.5], [0, n_genes - 1], color="#333", ls="--", lw=0.8, zorder=5)
    if tf_end > neural_end and tf_end < len(go_ids_filtered):
        ax.plot([tf_end - 0.5, tf_end - 0.5], [0, n_genes - 1], color="#333", ls="--", lw=0.8, zorder=5)

    ax.set_yticks(range(len(gene_ids)))
    ylabels = []
    for g in gene_ids:
        nm = label(neural, g)
        tag = "A" if proof_map.get(g, "") == "known_rnai_validated" else "B"
        sc = score_map.get(g, 0)
        ylabels.append(f"[{tag}] {nm} ({sc:.2f})")
    ax.set_yticklabels(ylabels, fontsize=5)
    for i, g in enumerate(gene_ids):
        if proof_map.get(g, "") == "known_rnai_validated":
            ax.get_yticklabels()[i].set_color(C_A)
            ax.get_yticklabels()[i].set_fontweight("bold")
    ax.set_ylabel("TF candidate ([A] Track A RNAi-validated, [B] Track B novel)", fontsize=7.5)

    ax.set_xticks(range(len(go_ids_filtered)))
    xlabels = [f"{term_name.get(t, t)[:24]} ({term_count.get(t, 0)})" for t in go_ids_filtered]
    ax.set_xticklabels(xlabels, rotation=60, ha="right", fontsize=5.5)
    for j, t in enumerate(go_ids_filtered):
        if term_neural.get(t, False):
            ax.get_xticklabels()[j].set_color(C_A)
            ax.get_xticklabels()[j].set_fontweight("bold")
        elif term_tf.get(t, False):
            ax.get_xticklabels()[j].set_color(C_B)

    legend_handles = [
        mpatches.Patch(facecolor=C_A, label="Annotated"),
        mpatches.Patch(facecolor="#F7F7F7", edgecolor="#CCCCCC", lw=0.5, label="Not annotated"),
        plt.Line2D([0], [0], color="#333", ls="--", lw=0.8, label="Group boundary"),
        mpatches.Patch(facecolor=C_A, alpha=0.3, label="Neural GO (+0.03 bonus)"),
        mpatches.Patch(facecolor=C_B, alpha=0.3, label="TF GO (+0.02 bonus)"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=6.5, frameon=True)
    ax.set_title("Gene Ontology (GO) annotation landscape across prioritized TF candidates",
                 fontsize=8.5, fontweight="bold", pad=12)

    fig.tight_layout()
    save_sup(fig, "fig_s1_go_landscape")
    print("  wrote fig_s1_go_landscape (PNG + PDF)")


def fig_s2_go_namespace_and_track():
    ref, mat = load_go()
    neural = load_neural()
    go_ids = [c for c in mat.columns if c.startswith("GO:")]
    term_neural = dict(zip(ref["go_id"], ref["neural_go"] == "yes"))
    term_ns = dict(zip(ref["go_id"], ref["namespace"]))
    term_count = dict(zip(ref["go_id"], ref["n_of_97_neural_candidates"]))
    go_meta = [(t, term_neural.get(t, False), term_ns.get(t, ""), term_count.get(t, 0))
               for t in go_ids if term_count.get(t, 0) >= 2]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_2COL, 2.8))

    # Panel a: Namespace distribution
    neural_ns, non_neural_ns = {}, {}
    for t, is_neural, ns, cnt in go_meta:
        ns_short = ns.split(" ")[0][:3] if ns else "unk"
        if is_neural:
            neural_ns[ns_short] = neural_ns.get(ns_short, 0) + 1
        else:
            non_neural_ns[ns_short] = non_neural_ns.get(ns_short, 0) + 1
    all_ns = sorted(set(list(neural_ns.keys()) + list(non_neural_ns.keys())))
    x = np.arange(len(all_ns))
    w = 0.35
    ax1.bar(x - w/2, [neural_ns.get(n, 0) for n in all_ns], w, color=C_A, alpha=0.85, label="Neural GO")
    ax1.bar(x + w/2, [non_neural_ns.get(n, 0) for n in all_ns], w, color=C_B, alpha=0.75, label="Non-neural GO")
    ax1.set_xticks(x)
    ns_labels = {"Bio": "Biological\nProcess", "Mol": "Molecular\nFunction",
                 "Cel": "Cellular\nComponent", "unk": "Unknown"}
    ax1.set_xticklabels([ns_labels.get(n, n) for n in all_ns], fontsize=7)
    ax1.set_ylabel("Number of GO terms", fontsize=7.5)
    ax1.legend(fontsize=6.5, frameon=False)
    panel_tag(ax1, "a")

    # Panel b: Track A vs Track B GO density
    proof_map = dict(zip(neural["gene_id"], neural["proof_status"]))
    track_a = [g for g in mat["gene_id"] if proof_map.get(g) == "known_rnai_validated"]
    track_b = [g for g in mat["gene_id"] if proof_map.get(g) == "novel_candidate"]
    mat_indexed = mat.set_index("gene_id")
    def per_track(gs):
        nn, oo = [], []
        for g in gs:
            if g in mat_indexed.index:
                row = mat_indexed.loc[g]
                nn.append(sum(1 for t in go_ids if row.get(t, 0) == 1 and term_neural.get(t, False)))
                oo.append(sum(1 for t in go_ids if row.get(t, 0) == 1 and not term_neural.get(t, False)))
        return nn, oo
    a_nn, a_oo = per_track(track_a)
    b_nn, b_oo = per_track(track_b)
    cats = ["Neural GO", "Other GO"]
    a_vals = [np.mean(a_nn) if a_nn else 0, np.mean(a_oo) if a_oo else 0]
    b_vals = [np.mean(b_nn) if b_nn else 0, np.mean(b_oo) if b_oo else 0]
    x2 = np.arange(len(cats))
    ax2.bar(x2 - w/2, a_vals, w, color=C_A, alpha=0.85, label=f"Track A (n={len(track_a)})")
    ax2.bar(x2 + w/2, b_vals, w, color=C_B, alpha=0.75, label=f"Track B (n={len(track_b)})")
    ax2.set_xticks(x2)
    ax2.set_xticklabels(cats, fontsize=7)
    ax2.set_ylabel("Mean GO terms per candidate", fontsize=7.5)
    ax2.legend(fontsize=6.5, frameon=False)
    panel_tag(ax2, "b")

    fig.tight_layout()
    save_sup(fig, "fig_s2_go_namespace_and_track")
    print("  wrote fig_s2_go_namespace_and_track (PNG + PDF)")


if __name__ == "__main__":
    fig_s1_go_landscape()
    fig_s2_go_namespace_and_track()

