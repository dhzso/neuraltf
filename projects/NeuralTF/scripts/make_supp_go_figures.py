"""Simplified, presentation-friendly supplementary figures for the NeuralTF
PlanMine GO-term data - designed so a reader can infer the message in seconds.

Figures (into ``projects/NeuralTF/figures/supplementary/``, 300 dpi):
- ``fig_s1_go_gene_term_map.png`` - 97 neural TFs (rows, ordered by track then
  integrated score) vs the most informative GO terms (<= 20 columns). Blue =
  annotated, white = not. No dendrograms, legible labels, one color legend.
- ``fig_s2_go_top10_dotmatrix.png`` - the final dual-track Top-10 x its key GO
  terms only. The "headline" figure: 10 rows, bold labels (A1..B5).
- ``fig_s3_go_top_terms.png`` - (a) top-15 GO terms by number of genes,
  (b) GO annotations per gene colored by proof status.
- ``fig_s4_go_neural_focus.png`` - neural-flagged GO terms only (<= 12 terms,
  most informative), same ordering as S1; the +0.03 neural-GO bonus base.
- ``go_gene_term_matrix_reduced.csv`` - binary matrix for a supplementary table.
- ``go_term_reference.csv`` (figures root) - every GO id used by the figures:
  canonical name, namespace, neural/TF flags, composite bonus, figure membership.

Uses the pipeline classifiers (go_term_flags, match_dna_binding_family) so the
figures always agree with the composite bonuses.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

try:
    from bioforge.projects.neuraltf.planmine import (
        go_term_flags,
        match_dna_binding_family,
    )
except ImportError:  # pragma: no cover - standalone fallback
    def go_term_flags(name):
        if not name:
            return False, False
        n = name.strip().lower()
        neural_kw = (
            "neuron", "nervous system development", "brain", "neurogenesis",
            "synaptic", "axon", "dendrite", "glial", "neural", "sensory",
            "auditory", "ophthalm", "eye ", "visual", "head", "cns",
        )
        tf_kw = (
            "transcription factor activity", "dna binding", "dna-binding",
            "regulation of transcription", "nucleic acid",
        )
        return any(k in n for k in neural_kw), any(k in n for k in tf_kw)

    def match_dna_binding_family(short_name):
        return []

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RUN = REPO_ROOT / "projects" / "NeuralTF" / "runs" / "pipeline_run"
DEFAULT_PARQUET = REPO_ROOT / "datasets" / "processed" / "planmine_annotations.parquet"
DEFAULT_TOP10 = REPO_ROOT / "projects" / "NeuralTF" / "results" / "top10_neural_tfs_prioritized.csv"
DEFAULT_OUT = REPO_ROOT / "projects" / "NeuralTF" / "figures" / "supplementary"

C_BLUE = "#0072B2"
C_ORANGE = "#E69F00"
C_SKY = "#56B4E9"
C_GRAY = "#999999"
C_PRESENT = "#0072B2"

PROOF_COLORS = {"known_rnai_validated": C_ORANGE, "novel_candidate": C_SKY}
PROOF_LABELS = {"known_rnai_validated": "RNAi-validated", "novel_candidate": "Novel"}
NS_LABELS = ["Biological process", "Molecular function",
             "Cellular component", "No namespace"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--run", type=Path, default=DEFAULT_RUN)
    p.add_argument("--parquet", type=Path, default=DEFAULT_PARQUET)
    p.add_argument("--top-csv", type=Path, default=DEFAULT_TOP10)
    p.add_argument("--obo", type=Path,
                   default=REPO_ROOT / "datasets" / "raw" / "go.obo",
                   help="go.obo ontology file (canonical names, namespaces)")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    return p.parse_args()


# ---------------------------------------------------------------------------
# Data loading & helpers
# ---------------------------------------------------------------------------

def _short_label(gene_id: str) -> str:
    parts = gene_id.split("_")
    if gene_id.startswith("dd_Smed_") and len(parts) > 3:
        return "dd" + parts[3]
    return gene_id


def load_obo(path: str | Path) -> dict[str, tuple[str, str, bool]]:
    """Parse go.obo -> {id: (canonical_name, namespace, is_obsolete)}.

    Includes alt_id aliases so every PlanMine term can be resolved."""
    path = Path(path)
    if not path.exists():
        print(f"  [warn] go.obo not found at {path} - falling back to "
              "PlanMine's own term names. Download the current release from "
              "https://current.geneontology.org/ontology/go.obo and place it "
              "at datasets/raw/go.obo (see datasets/MANIFEST.md).")
        return {}
    out: dict[str, tuple[str, str, bool]] = {}
    cur_id, cur_name, cur_ns, cur_obs = None, "", "", False
    for line in Path(path).open(encoding="utf-8", errors="replace"):
        line = line.strip()
        if line == "[Term]":
            if cur_id is not None:
                out[cur_id] = (cur_name, cur_ns, cur_obs)
            cur_id, cur_name, cur_ns, cur_obs = None, "", "", False
        elif line.startswith("id: "):
            cur_id = line[4:]
        elif line.startswith("name: "):
            cur_name = line[6:]
        elif line.startswith("namespace: "):
            cur_ns = line[11:]
        elif line.startswith("is_obsolete: true"):
            cur_obs = True
        elif line.startswith("alt_id: ") and cur_id:
            out[line[8:]] = (cur_name, cur_ns, cur_obs)
    if cur_id is not None:
        out[cur_id] = (cur_name, cur_ns, cur_obs)
    return out


def load_data(args) -> tuple:
    """Return (go rows with flags, labels, names, neural df, top df, ids)."""
    ann = pd.read_parquet(args.parquet)
    neural = pd.read_csv(args.run / "rank_neural.csv")
    ids = sorted(neural["gene_id"].tolist())
    id_set = set(ids)

    go = ann[(ann["kind"] == "go") & (ann["gene_id_v6"].isin(id_set))].copy()
    go = go.drop_duplicates(["gene_id_v6", "key"])
    go["annotated"] = 1  # presence marker (the matrix must record annotation, not flags)

    obo = load_obo(args.obo)
    meta = {t: obo.get(t, ("", "", False)) for t in go["key"].unique()}
    go["obo_name"] = [meta[t][0] for t in go["key"]]
    go["obo_namespace"] = [meta[t][1] for t in go["key"]]
    go["is_obsolete"] = [meta[t][2] for t in go["key"]]

    fallback = go["value"].fillna("").astype(str)
    go["name"] = [meta[t][0] or f for t, f in zip(go["key"], fallback)]
    go["namespace"] = [meta[t][1] or n for t, n in
                       zip(go["key"], go["namespace"].fillna("").astype(str))]
    go["is_neural"], go["is_tf"] = zip(*go["name"].map(go_term_flags))

    labels = {}
    names = {}
    for gid in ids:
        labels[gid] = _short_label(gid)
        nm = neural.loc[neural["gene_id"] == gid, "gene_name"]
        val = str(nm.iloc[0]) if not nm.isna().all() else ""
        names[gid] = val if val not in ("nan", "") else labels[gid]

    top = pd.read_csv(args.top_csv) if args.top_csv.exists() else None
    if top is not None:
        for _, r in top.iterrows():
            nm = str(r["gene_name"])
            if nm not in ("nan", ""):
                names[r["gene_id_v6"]] = nm
    return go, labels, names, neural, top, ids


def order_genes(ids: list[str], neural: pd.DataFrame, top) -> list[str]:
    """Track A/B genes first (by composite), then the rest by integrated."""
    score = {g: s for g, s in zip(neural["gene_id"], neural["integrated_score"])}
    order: list[str] = []
    if top is not None:
        for rank_src in ("A", "B"):
            sub = top[top["track"] == rank_src].sort_values(
                "composite_score", ascending=False)
            order += sub["gene_id_v6"].tolist()
    rest = [g for g in ids if g not in order]
    rest.sort(key=lambda g: score.get(g, 0.0), reverse=True)
    return order + rest


def track_no(gene_id: str, top) -> str | None:
    if top is None:
        return None
    row = top[top["gene_id_v6"] == gene_id]
    if row.empty:
        return None
    return f"{row.iloc[0]['track']}{row.iloc[0]['rank']}"


def pick_terms(go: pd.DataFrame, n_neural: int = 8, n_top: int = 12,
               neural_only: bool = False, max_terms: int = 20) -> list[str]:
    """Most informative non-obsolete terms: neural-flagged first, then most shared."""
    cnt = go.groupby("key").size()
    neural_of = go.groupby("key")["is_neural"].any()
    obsolete = go.groupby("key")["is_obsolete"].any()
    keys = sorted(cnt.index.tolist(),
                  key=lambda t: (bool(obsolete[t]), not bool(neural_of[t]),
                                 -int(cnt[t]), t))
    if neural_only:
        keys = [t for t in keys if neural_of[t] and not obsolete[t]]
        return keys[:max_terms]
    keep = [t for t in keys if neural_of[t]][:n_neural]
    keep += [t for t in keys if not neural_of[t]][:n_top]
    return keep[:max_terms]


def build_matrix(go: pd.DataFrame, ids: list[str], terms: list[str],
                 gene_order: list[str]) -> pd.DataFrame:
    """Binary gene x term matrix; cell = 1 iff the gene is annotated with the term."""
    mat = go.pivot_table(index="gene_id_v6", columns="key", values="annotated",
                         aggfunc="first").reindex(gene_order).fillna(0).astype(int)
    return mat[[t for t in terms if t in mat.columns]]


def term_display(go: pd.DataFrame, terms: list[str]) -> dict[str, str]:
    """Short display names per term; duplicate names get a '*' marker."""
    name_of = go.groupby("key")["name"].apply(
        lambda s: s.mode().iloc[0] if not s.mode().empty else "")
    out = {}
    seen = set()
    for t in terms:
        nm = str(name_of.get(t, t)) or t
        if nm in seen:
            nm = nm + " *"
        seen.add(nm)
        out[t] = nm
    return out


def _draw_map(ax, mat: pd.DataFrame, display_terms: dict[str, str], gene_labels,
              row_strip_colors, strip_label, fig_height: float = 12.5) -> None:
    """Simple binary map: one colored square per annotation, white grid lines."""
    import numpy as np
    from matplotlib.colors import ListedColormap
    data = mat.to_numpy()
    n_rows, n_cols = data.shape
    ax.imshow(data, cmap=ListedColormap(["#ffffff", C_PRESENT]), vmin=0, vmax=1,
              aspect="auto", interpolation="nearest")
    for i in range(n_rows + 1):
        ax.axhline(i - 0.5, color="white", lw=1.2)
    for j in range(n_cols + 1):
        ax.axvline(j - 0.5, color="white", lw=1.2)

    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(gene_labels, fontsize=8)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([display_terms[t] for t in mat.columns],
                       fontsize=6.5, rotation=45, ha="right")
    ax.tick_params(left=False, bottom=False)

    if row_strip_colors is not None:
        xs = np.full(n_rows, n_cols + 0.55)
        ys = np.arange(n_rows)
        s = (fig_height * 72.0 / max(n_rows, 1)) ** 2 * 0.9
        ax.scatter(xs, ys, c=row_strip_colors, marker="s", s=s,
                   edgecolors="none", clip_on=False, zorder=5)
        ax.text(n_cols + 0.95, -1.4, strip_label, fontsize=7, color="black",
                ha="left", va="bottom", clip_on=False)


def _save(fig, out: Path, name: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.png", facecolor="white", bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  wrote {name}.png")


def _legend(ax, handles, loc: str = "lower left", ncol: int = 1) -> None:
    ax.legend(handles=handles, loc=loc, fontsize=8, frameon=False)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def fig_s1_map(go, neural, ids, labels, names, top, out: Path) -> None:
    terms = pick_terms(go, n_neural=8, n_top=12)
    gene_order = order_genes(ids, neural, top)
    mat = build_matrix(go, gene_order, terms, gene_order)
    display = term_display(go, mat.columns)
    display_map = dict(zip(mat.columns, display))

    proof = {g: p for g, p in zip(neural["gene_id"], neural["proof_status"])}
    strip = [PROOF_COLORS.get(str(proof.get(g, "")), C_GRAY) for g in gene_order]
    score = {g: s for g, s in zip(neural["gene_id"], neural["integrated_score"])}

    fig, ax = plt.subplots(figsize=(11, 12.5))
    labels_ = []
    for i, g in enumerate(gene_order):
        tn = track_no(g, top)
        pre = f"{tn}·" if tn else ""
        sc = f" ({score.get(g, float('nan')):.3f})" if i < 10 else ""
        labels_.append(f"{pre}{names[g]}{sc}")

    _draw_map(ax, mat, display_map, labels_, strip, "proof status")
    ax.set_xlim(-0.5, mat.shape[1] + 2.2)
    ax.set_title(
        f"PlanMine GO annotations of the {len(gene_order)} neural TF candidates\n"
        f"Blue = the gene is annotated with that GO term (top informative terms shown); "
        f"rows ordered by track then score (top-10 rows show integrated score)",
        fontsize=10, fontweight="bold", pad=20)
    handles = [mpatches.Patch(color=PROOF_COLORS[k], label=v)
                     for k, v in PROOF_LABELS.items()]
    _legend(ax, handles, loc="lower left")
    _save(fig, out, "fig_s1_go_gene_term_map")


def fig_s2_top10_dotmatrix(go, neural, ids, labels, names, top, out: Path) -> None:
    if top is None:
        return
    terms = pick_terms(go, n_neural=6, n_top=9, max_terms=15)
    gene_order = [g for g in order_genes(ids, neural, top) if g in set(top["gene_id_v6"])]
    mat = build_matrix(go, gene_order, terms, gene_order)
    display = term_display(go, mat.columns)
    display_map = dict(zip(mat.columns, display))

    fig, ax = plt.subplots(figsize=(10.5, 6))
    labels_ = []
    for g in gene_order:
        tn = track_no(g, top)
        labels_.append(f"{tn}  {names[g]}")
    _draw_map(ax, mat, display_map, labels_, None, None)
    for j in range(mat.shape[1]):
        ax.add_patch(plt.Rectangle((j - 0.5, mat.shape[0] - 0.5), 1.0, 0.5,
                                   facecolor=C_GRAY, alpha=0.25,
                                   edgecolor="none", clip_on=False))
    ax.set_title(
        "GO-term profile of the final dual-track Top-10\n"
        "blue = PlanMine GO annotation present",
        fontsize=10.5, fontweight="bold", pad=14)
    handles = [mpatches.Patch(color=C_PRESENT, label="GO annotation present"),
               mpatches.Patch(color=C_GRAY, alpha=0.25, label="no annotation")]
    _legend(ax, handles, loc="lower left")
    _save(fig, out, "fig_s2_go_top10_dotmatrix")


def fig_s3_top_terms(go, neural, ids, labels, names, top, out: Path) -> None:
    cnt = go.groupby("key").size().sort_values(ascending=False)
    neural_of = go.groupby("key")["is_neural"].any()

    top_terms = cnt.head(15).iloc[::-1]
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.2),
                             gridspec_kw={"width_ratios": [1.7, 1.0]})
    ax = axes[0]
    colors = [C_ORANGE if neural_of[t] else "#7F7F7F" for t in top_terms.index]
    ax.barh(range(len(top_terms)), top_terms.values, color=colors, height=0.72)
    ax.set_yticks(range(len(top_terms)))
    ax.set_yticklabels([t.replace("GO:", "") for t in top_terms.index],
                       fontsize=7.5)
    for i, v in enumerate(top_terms.values):
        ax.text(v + 0.15, i, str(int(v)), va="center", fontsize=7.5)
    ax.set_xlabel("Number of the 97 neural TF candidates")
    ax.set_xlim(0, max(top_terms) + 3)
    ax.tick_params(direction="out")
    ax.set_title("Top GO terms (orange = neural-related)", fontsize=9,
                 fontweight="bold", pad=8)

    ax2 = axes[1]
    per_gene = go.groupby("gene_id_v6").size()
    proof = {g: p for g, p in zip(neural["gene_id"], neural["proof_status"])}
    order = per_gene.sort_values(ascending=False)
    ax2.bar(range(len(order)), order.values,
            color=[PROOF_COLORS.get(str(proof.get(g, "")), C_GRAY)
                   for g in order.index], width=0.8)
    n_no_go = sum(1 for g in ids if g not in per_gene.index)
    ax2.set_xlabel("Genes (sorted by annotation count)")
    ax2.set_ylabel("GO annotations")
    ax2.set_ylim(0, order.max() + 9)
    ax2.set_xticks([])
    ax2.tick_params(direction="out")
    stats = (f"{len(order)}/{len(ids)} genes have \u22651 GO term; "
             f"median {int(per_gene.median())} per gene")
    ax2.text(0.02, 0.97, stats, transform=ax2.transAxes, ha="left", va="top",
             fontsize=7.5, color=C_GRAY)
    ax2.set_title("GO annotations per candidate", fontsize=9,
                  fontweight="bold", pad=8)
    handles = [mpatches.Patch(color=PROOF_COLORS[k], label=v)
                     for k, v in PROOF_LABELS.items()]
    _legend(axes[1], handles, loc="upper right")
    _save(fig, out, "fig_s3_go_top_terms")


def fig_s4_neural_focus(go, neural, ids, labels, names, top, out: Path) -> None:
    terms = pick_terms(go, neural_only=True, max_terms=12)
    if len(terms) < 2:
        print("  fig_s4: fewer than 2 neural terms, skipped")
        return
    gene_order = order_genes(ids, neural, top)
    mat = build_matrix(go, gene_order, terms, gene_order)
    display = term_display(go, mat.columns)
    display_map = dict(zip(mat.columns, display))

    proof = {g: p for g, p in zip(neural["gene_id"], neural["proof_status"])}
    strip = [PROOF_COLORS.get(str(proof.get(g, "")), C_GRAY) for g in gene_order]

    fig, ax = plt.subplots(figsize=(9.5, 12.5))
    labels_ = []
    for g in gene_order:
        tn = track_no(g, top)
        labels_.append(f"{tn}·" if tn else "")
    labels_ = [f"{pre}{names[g]}" for pre, g in zip(labels_, gene_order)]

    _draw_map(ax, mat, display_map, labels_, strip, "proof status")
    ax.set_xlim(-0.5, mat.shape[1] + 2.2)
    ax.set_title(
        f"Neural-related GO terms only ({mat.shape[1]} terms) -\n"
        f"the evidence behind the +0.03 neural-GO composite bonus",
        fontsize=10, fontweight="bold", pad=20)
    handles = [mpatches.Patch(color=PROOF_COLORS[k], label=v)
                     for k, v in PROOF_LABELS.items()]
    _legend(ax, handles, loc="lower left")
    _save(fig, out, "fig_s4_go_neural_focus")


def write_matrix_csv(go, ids, top, out: Path) -> None:
    terms = pick_terms(go)
    mat = build_matrix(go, ids, terms, ids).copy()
    cell = []
    for g in mat.index:
        nm = ""
        if top is not None:
            row = top[top["gene_id_v6"] == g]
            nm = str(row.iloc[0]["gene_name"]) if not row.empty else ""
        cell.append(nm if nm not in ("nan", "") else _short_label(g))
    mat.insert(0, "cell", cell)
    mat.to_csv(out / "go_gene_term_matrix_reduced.csv", index_label="gene_id")
    print(f"  wrote go_gene_term_matrix_reduced.csv ({mat.shape[0]} x {mat.shape[1]})")


def write_go_term_reference(go, ids, top, out: Path) -> None:
    """Write figures/go_term_reference.csv — what each GO id used by the
    figures corresponds to: canonical name, namespace, neural/TF flags, the
    composite bonus it triggers, and which figure shows the term."""
    cnt = go.groupby("key").size()
    neural_of = go.groupby("key")["is_neural"].any()
    tf_of = go.groupby("key")["is_tf"].any()
    obsolete_of = go.groupby("key")["is_obsolete"].any()
    name_of = go.groupby("key")["name"].apply(
        lambda s: s.mode().iloc[0] if not s.mode().empty else "")
    ns_of = go.groupby("key")["namespace"].apply(
        lambda s: s.mode().iloc[0] if not s.mode().empty else "")

    s1 = set(pick_terms(go, n_neural=8, n_top=12))
    s2 = set(pick_terms(go, n_neural=6, n_top=9, max_terms=15))
    s3 = set(cnt.sort_values(ascending=False).head(15).index)
    s4 = set(pick_terms(go, neural_only=True, max_terms=12))
    s9: set[str] = set()
    if top is not None:
        top_ids = top["gene_id_v6"].astype(str).tolist()
        sub = go[go["gene_id_v6"].isin(top_ids)].groupby("key").size()
        s9 = set(sub[sub >= 2].sort_values(ascending=False).head(12).index)

    rows = []
    for t in sorted(cnt.index,
                    key=lambda k: (not bool(neural_of[k]), -int(cnt[k]), k)):
        bonus = []
        if neural_of[t]:
            bonus.append("go_neural +0.03")
        if tf_of[t]:
            bonus.append("go_tf +0.02")
        rows.append({
            "go_id": t,
            "term": name_of.get(t, ""),
            "namespace": ns_of.get(t, ""),
            "neural_go": "yes" if neural_of[t] else "no",
            "tf_go": "yes" if tf_of[t] else "no",
            "composite_bonus": " + ".join(bonus),
            "n_of_97_neural_candidates": int(cnt[t]),
            "in_fig9_go_dotplot": "yes" if t in s9 else "no",
            "in_fig_s1_map": "yes" if t in s1 else "no",
            "in_fig_s2_top10_matrix": "yes" if t in s2 else "no",
            "in_fig_s3_top15": "yes" if t in s3 else "no",
            "in_fig_s4_neural_focus": "yes" if t in s4 else "no",
            "is_obsolete": "yes" if obsolete_of[t] else "no",
        })
    ref = pd.DataFrame(rows)
    ref.to_csv(out / "go_term_reference.csv", index=False)
    print(f"  wrote go_term_reference.csv ({len(ref)} terms)")


def main() -> None:
    args = parse_args()
    go, labels, names, neural, top, ids = load_data(args)
    out = args.out
    out.mkdir(parents=True, exist_ok=True)

    print("Figures:")
    fig_s1_map(go, neural, ids, labels, names, top, out)
    fig_s2_top10_dotmatrix(go, neural, ids, labels, names, top, out)
    fig_s3_top_terms(go, neural, ids, labels, names, top, out)
    fig_s4_neural_focus(go, neural, ids, labels, names, top, out)
    write_matrix_csv(go, ids, top, out)
    write_go_term_reference(go, ids, top, out.parent)


if __name__ == "__main__":
    main()