"""GO term figures — all candidates x GO terms with neural highlighted, top-10 profiles, namespace analysis."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, matplotlib.patches as mpatches
import numpy as np, pandas as pd

SUP = FIG / "supplementary"
SUP.mkdir(parents=True, exist_ok=True)
GO_REF = FIG / "go_term_reference.csv"
GO_MATRIX = SUP / "go_gene_term_matrix_reduced.csv"
C_NEURAL_GO = "#D55E00"
C_TF_GO = "#0072B2"
C_OTHER_GO = "#999999"

def load_go():
    return pd.read_csv(GO_REF), pd.read_csv(GO_MATRIX)


def fig_s5_go_heatmap_all99():
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

    # Find boundaries between groups
    neural_end = sum(1 for t in go_ids_filtered if term_neural.get(t, False))
    tf_end = neural_end + sum(1 for t in go_ids_filtered[neural_end:]
                              if term_tf.get(t, False))

    proof_map = dict(zip(neural["gene_id"], neural["proof_status"]))
    score_map = dict(zip(neural["gene_id"], neural["integrated_score"]))
    gene_ids = sorted(mat["gene_id"].tolist(),
                      key=lambda g: (proof_map.get(g, "") != "known_rnai_validated",
                                     -score_map.get(g, 0)))
    mat_indexed = mat.set_index("gene_id").loc[gene_ids]
    data = mat_indexed[go_ids_filtered].fillna(0).values

    fig, ax = plt.subplots(figsize=(14, 13))

    # Professional color scheme: soft coral for annotated, light gray for absent
    cmap = plt.cm.colors.ListedColormap(["#F5F5F5", "#C44E52"])
    ax.imshow(data, cmap=cmap, aspect="auto", interpolation="nearest",
              vmin=0, vmax=1)

    # Dotted vertical lines separating groups — only within heatmap area
    n_genes = len(gene_ids)
    if neural_end > 0 and neural_end < len(go_ids_filtered):
        ax.plot([neural_end - 0.5, neural_end - 0.5], [0, n_genes - 1],
                color="#333", ls="--", lw=1.2, zorder=5)
    if tf_end > neural_end and tf_end < len(go_ids_filtered):
        ax.plot([tf_end - 0.5, tf_end - 0.5], [0, n_genes - 1],
                color="#333", ls="--", lw=1.2, zorder=5)

    # Group labels below x-tick labels
    if neural_end > 0:
        ax.text(neural_end / 2 - 0.5, -3.0, "Neural GO", fontsize=9,
                ha="center", fontweight="bold", color="#C44E52")
    if tf_end > neural_end:
        ax.text((neural_end + tf_end) / 2 - 0.5, -3.0, "TF GO", fontsize=9,
                ha="center", fontweight="bold", color="#4C72B0")
    remaining = len(go_ids_filtered) - tf_end
    if remaining > 0:
        ax.text(tf_end + remaining / 2 - 0.5, -3.0, "Other GO", fontsize=9,
                ha="center", fontweight="bold", color="#555555")

    # Y-axis labels with decoded information
    ax.set_yticks(range(len(gene_ids)))
    ylabels = []
    for g in gene_ids:
        nm = label(neural, g)
        tag = "A" if proof_map.get(g, "") == "known_rnai_validated" else "B"
        sc = score_map.get(g, 0)
        ylabels.append(f"[{tag}] {nm} ({sc:.3f})")
    ax.set_yticklabels(ylabels, fontsize=5.5)
    for i, g in enumerate(gene_ids):
        if proof_map.get(g, "") == "known_rnai_validated":
            ax.get_yticklabels()[i].set_color("#C44E52")
    ax.set_ylabel("TF candidate  —  [A] Track A (RNAi-validated),  [B] Track B (novel)\n"
                   "Integrated evidence score in parentheses",
                   fontsize=8, fontweight="bold", labelpad=12)

    # X-axis labels
    ax.set_xticks(range(len(go_ids_filtered)))
    xlabels = [f"{term_name.get(t, t)}\n({term_count.get(t, 0)})" for t in go_ids_filtered]
    ax.set_xticklabels(xlabels, rotation=60, ha="right", fontsize=5.5)
    for j, t in enumerate(go_ids_filtered):
        if term_neural.get(t, False):
            ax.get_xticklabels()[j].set_color("#C44E52")
            ax.get_xticklabels()[j].set_fontweight("bold")
        elif term_tf.get(t, False):
            ax.get_xticklabels()[j].set_color("#4C72B0")
        else:
            ax.get_xticklabels()[j].set_color("#555555")
    ax.set_xlabel("GO term name  (number of candidates annotated with this term)",
                   fontsize=8, fontweight="bold", labelpad=40)

    # Legend
    legend_handles = [
        mpatches.Patch(facecolor="#C44E52", label="Annotated (present)"),
        mpatches.Patch(facecolor="#F5F5F5", edgecolor="#CCC", label="Not annotated"),
        plt.Line2D([0],[0], color="#333", ls="--", lw=1.2, label="Group boundary"),
        mpatches.Patch(facecolor="#C44E52", alpha=0.3, label="Neural GO (+0.03 bonus)"),
        mpatches.Patch(facecolor="#4C72B0", alpha=0.3, label="TF GO (+0.02 bonus)"),
        mpatches.Patch(facecolor="#C44E52", label="Track A (RNAi-validated)"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=7.5, frameon=True,
              title="Legend", title_fontsize=9)

    ax.set_title("GO annotation landscape — which GO terms are shared across neural TF candidates?\n"
                 "Neural GO terms linked to brain development, TF GO terms linked to transcription regulation, "
                 "Other GO terms = remaining functional annotations",
                 fontweight="bold", pad=15, fontsize=9)
    ax.set_ylim(len(gene_ids) - 0.5, -3.5)
    ax.set_xlim(-0.5, len(go_ids_filtered) - 0.5)

    fig.tight_layout()
    fig.savefig(SUP / "fig_s5_go_heatmap_neural.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  wrote fig_s5_go_heatmap_neural.png")
    plt.close(fig)
    print("  wrote fig_s5_go_heatmap_neural.png")


def fig_s6_top10_go_profiles():
    ref, mat = load_go()
    top10 = load_top10()
    neural = load_neural()
    go_ids = [c for c in mat.columns if c.startswith("GO:")]
    term_neural = dict(zip(ref["go_id"], ref["neural_go"] == "yes"))
    term_tf = dict(zip(ref["go_id"], ref["tf_go"] == "yes"))
    term_count = dict(zip(ref["go_id"], ref["n_of_97_neural_candidates"]))
    go_ids_filtered = [t for t in go_ids if term_count.get(t, 0) >= 2]
    mat_indexed = mat.set_index("gene_id")
    n = len(top10)
    fig, axes = plt.subplots(2, 5, figsize=(18, 8), sharey=False)
    axes = axes.flatten()
    for idx, (_, row) in enumerate(top10.iterrows()):
        gid = row["gene_id"]
        ax = axes[idx]
        nm = label(neural, gid)
        track = row.get("track", "")
        if gid in mat_indexed.index:
            gene_go = mat_indexed.loc[gid]
            annotations = [t for t in go_ids_filtered if gene_go.get(t, 0) == 1]
        else:
            annotations = []
        neural_count = sum(1 for t in annotations if term_neural.get(t, False))
        tf_count = sum(1 for t in annotations if term_tf.get(t, False) and not term_neural.get(t, False))
        other_count = len(annotations) - neural_count - tf_count
        no_go = len(go_ids_filtered) - len(annotations)
        counts = [neural_count, tf_count, other_count, no_go]
        labels_ = ["Neural GO", "TF GO", "Other GO", "Not annotated"]
        colors = ["#C44E52", "#4C72B0", "#55A868", "#E0E0E0"]
        ax.barh(range(4), counts, color=colors, edgecolor="white", lw=0.5, height=0.6)
        ax.set_yticks(range(4))
        ax.set_yticklabels(labels_, fontsize=6)
        for i, v in enumerate(counts):
            if v > 0:
                ax.text(v + 0.1, i, str(v), fontsize=7, va="center", fontweight="bold")
        tc = "#C44E52" if track == "A" else "#4C72B0"
        ax.set_title(f"[{track}] {nm}", fontsize=8, fontweight="bold", color=tc, pad=5)
        ax.set_xlabel("GO terms", fontsize=6)
        ax.set_xlim(0, max(counts) + 2)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.suptitle("GO annotation profiles — Top 10 candidates\n"
                 "Red = neural GO, Blue = TF GO, Green = other GO",
                 fontweight="bold", fontsize=10, y=1.02)
    fig.tight_layout(w_pad=1)
    fig.savefig(SUP / "fig_s6_top10_go_profiles.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  wrote fig_s6_top10_go_profiles.png")


def fig_s7_go_namespace_and_track():
    ref, mat = load_go()
    neural = load_neural()
    go_ids = [c for c in mat.columns if c.startswith("GO:")]
    term_neural = dict(zip(ref["go_id"], ref["neural_go"] == "yes"))
    term_ns = dict(zip(ref["go_id"], ref["namespace"]))
    term_count = dict(zip(ref["go_id"], ref["n_of_97_neural_candidates"]))
    go_meta = [(t, term_neural.get(t, False), term_ns.get(t, ""), term_count.get(t, 0))
               for t in go_ids if term_count.get(t, 0) >= 2]
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax = axes[0]
    neural_ns, non_neural_ns = {}, {}
    for t, is_neural, ns, cnt in go_meta:
        ns_short = ns.split(" ")[0][:3] if ns else "unk"
        if is_neural:
            neural_ns[ns_short] = neural_ns.get(ns_short, 0) + 1
        else:
            non_neural_ns[ns_short] = non_neural_ns.get(ns_short, 0) + 1
    all_ns = sorted(set(list(neural_ns.keys()) + list(non_neural_ns.keys())))
    x = np.arange(len(all_ns)); w = 0.35
    ax.bar(x - w/2, [neural_ns.get(n, 0) for n in all_ns], w, color="#C44E52", alpha=0.8, label="Neural GO")
    ax.bar(x + w/2, [non_neural_ns.get(n, 0) for n in all_ns], w, color="#4C72B0", alpha=0.6, label="Non-neural GO")
    ax.set_xticks(x)
    ns_labels = {"Bio": "Biological\nProcess", "Mol": "Molecular\nFunction",
                 "Cel": "Cellular\nComponent", "unk": "Unknown"}
    ax.set_xticklabels([ns_labels.get(n, n) for n in all_ns], fontsize=8)
    ax.set_ylabel("Number of GO terms"); ax.set_title("GO namespace distribution", fontweight="bold", pad=8)
    ax.legend(fontsize=7, frameon=True); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax2 = axes[1]
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
    x = np.arange(len(cats)); w = 0.35
    ax2.bar(x - w/2, a_vals, w, color="#C44E52", alpha=0.8, label=f"Track A (n={len(track_a)})")
    ax2.bar(x + w/2, b_vals, w, color="#4C72B0", alpha=0.8, label=f"Track B (n={len(track_b)})")
    ax2.set_xticks(x); ax2.set_xticklabels(cats, fontsize=8)
    ax2.set_ylabel("Mean GO annotations per candidate")
    ax2.set_title("GO density by track", fontweight="bold", pad=8)
    ax2.legend(fontsize=7, frameon=True); ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)
    fig.suptitle("GO term analysis — namespace and track comparison", fontweight="bold", fontsize=11, y=1.02)
    fig.tight_layout(w_pad=2)
    fig.savefig(SUP / "fig_s7_go_namespace_track.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  wrote fig_s7_go_namespace_track.png")


if __name__ == "__main__":
    fig_s5_go_heatmap_all99()
    fig_s6_top10_go_profiles()
    fig_s7_go_namespace_and_track()
