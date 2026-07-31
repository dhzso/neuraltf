"""Built-in workflow steps wiring Layers 8A / 8B to the Layer 7 registry.

This module is imported lazily by the CLI / UI so the framework's omics
dependencies don't have to be importable when only the registry is used.
Every callable here has a small docstring so an LLM-driven assistant can
introspect workflows through the StepRegistry.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from bioforge.core.logging import get_logger
from bioforge.evidence.schema import EvidenceRecord, EvidenceSource
from bioforge.ingest import ingest_dataset
# Lazy import — batch step imports run_harmony only when called
from bioforge.omics.cluster import run_cluster
from bioforge.omics.normalize import run_normalize
from bioforge.omics.qc import run_qc
from bioforge.workflow.registry import register

logger = get_logger("workflow.steps")


@register("ingest", outputs=["adata"],
          description="Auto-detect dataset format and load as AnnData.")
def ingest(source: str, **kwargs) -> dict:
    adata = ingest_dataset(source, **kwargs)
    return {"adata": adata}


@register("qc", inputs=["adata"], outputs=["adata"],
          description="Run QC pipeline (filter + metrics) on an AnnData.")
def qc(adata: ad.AnnData, **kwargs) -> dict:
    return {"adata": run_qc(adata, **kwargs)}


@register("normalize", inputs=["adata"], outputs=["adata"],
          description="Normalize + log1p + HVG.")
def normalize(adata: ad.AnnData, **kwargs) -> dict:
    return {"adata": run_normalize(adata, **kwargs)}


@register("cluster", inputs=["adata"], outputs=["adata"],
          description="PCA + neighbors + leiden + UMAP.")
def cluster(adata: ad.AnnData, **kwargs) -> dict:
    return {"adata": run_cluster(adata, **kwargs)}


@register("trajectory", inputs=["adata"], outputs=["adata"],
          description="PAGA trajectory inference (cluster graph abstraction).")
def trajectory(adata: ad.AnnData, **kwargs) -> dict:
    from bioforge.omics.trajectory import paga
    return {"adata": paga(adata, **kwargs)}


@register("batch_correct", inputs=["adata"], outputs=["adata"],
          description="Harmony batch correction.")
def batch_correct(adata: ad.AnnData, **kwargs) -> dict:
    from bioforge.omics.batch import run_harmony
    return {"adata": run_harmony(adata, **kwargs)}


@register(
    "evidence.demo_rank",
    outputs=["records"],
    description="Build a tiny demo EvidenceRecord list so workflows can run \
                 end-to-end without external datasets.",
)
def demo_rank() -> dict:
    r1 = EvidenceRecord(gene_id="dd_Smed_v6_x1", gene_name="soxB")
    r1.add_score(EvidenceSource.EXPRESSION, 1.0)
    r1.add_score(EvidenceSource.SPECIFICITY, 0.8)
    r1.add_score(EvidenceSource.REPRODUCIBILITY, 1.0)
    r1.add_score(EvidenceSource.RNai, 0.0)
    r1.add_score(EvidenceSource.CORRELATION, 0.6)
    r2 = EvidenceRecord(gene_id="dd_Smed_v6_x2", gene_name="newCandidate")
    r2.add_score(EvidenceSource.EXPRESSION, 0.9)
    r2.add_score(EvidenceSource.SPECIFICITY, 0.95)
    r2.add_score(EvidenceSource.REPRODUCIBILITY, 0.66)
    r2.add_score(EvidenceSource.RNai, 0.0)
    records = [r1, r2]
    return {"records": records}


@register(
    "evidence.write_rank_csv",
    inputs=["records"],
    outputs=["csv_path"],
    description="Write EvidenceRecord list as CSV.",
)
def write_rank_csv(records: list[EvidenceRecord], path: str = "rank.csv") -> dict:
    rows = []
    from bioforge.evidence import EvidenceScorer
    s = EvidenceScorer()
    for r in records:
        rows.append({
            "gene_id": r.gene_id,
            "gene_name": r.gene_name or "",
            "integrated_score": s.integrated_score(r),
            "n_streams": r.supporting_streams(),
            **{src.value: r.scores.get(src, None) for src in EvidenceSource},
        })
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    return {"csv_path": path}


@register(
    "evidence.build_cards",
    inputs=["records"],
    outputs=["cards"],
    description="Build EvidenceCard objects from EvidenceRecord list.",
)
def build_cards(records: list[EvidenceRecord], **kwargs) -> dict:
    from bioforge.evidence import build_cards_for_records
    atlas_membership = kwargs.pop("atlas_membership", None)
    prior_fstf_ids = kwargs.pop("prior_fstf_ids", None)
    if isinstance(atlas_membership, dict):
        atlas_membership = {k: set(v) for k, v in atlas_membership.items()}
    if isinstance(prior_fstf_ids, list):
        prior_fstf_ids = set(prior_fstf_ids)
    cards = build_cards_for_records(
        records,
        atlas_membership=atlas_membership,
        prior_fstf_ids=prior_fstf_ids,
    )
    return {"cards": cards}


@register(
    "report.write_cards_md",
    inputs=["cards"],
    outputs=["md_path"],
    description="Render EvidenceCards as a single markdown report.",
)
def write_cards_md(cards: list, path: str = "evidence_cards.md") -> dict:
    from bioforge.evidence.cards import render_cards_markdown
    md = render_cards_markdown(cards)
    Path(path).write_text(md, encoding="utf-8")
    return {"md_path": path}


@register(
    "ai.summarize_candidates",
    inputs=["cards"],
    outputs=["summary"],
    description="Use the AI assistant (or StubAssistant fallback) to summarize candidate TFs.",
)
def summarize_candidates(cards: list, **kwargs) -> dict:
    from bioforge.ai import build_assistant
    from bioforge.ai.assistant import ChatMessage
    assistant = build_assistant()
    # Compact representation for the assistant
    summary_lines = [
        f"- {c.gene_name or c.gene_id}: tier={c.tier.value}, "
        f"score={c.integrated_score:.2f}, status={c.proof_status.value}, "
        f"atlases={sorted(c.atlases_supported)}"
        for c in cards
    ]
    prompt = (
        "I have ranked these planarian TF candidates by evidence. Write a "
        "short paragraph naming the top novel (not-yet-RNAi-tested) TFs to "
        "prioritate, separate from any already RNAi-validated, and note what "
        "I should chase next experimentally.\n\n" + "\n".join(summary_lines)
    )
    resp = assistant.complete([ChatMessage(role="user", content=prompt)])
    return {"summary": resp.content, "model": resp.model}


# ---------------------------------------------------------------------------
# Real-data NeuralTF steps — registered under ingest.fincher / ingest.plass /
# ingest.king and evidence.* targets so projects can express the full pipeline
# in YAML. Reader targets are dataset-specific but return AnnData / DataFrames
# so the loop construct can iterate over them uniformly.
# ---------------------------------------------------------------------------


@register("ingest.fincher", outputs=["adata"],
          description="Load the Fincher 2018 principal DGE matrix as AnnData (dd_Smed_v4).")
def ingest_fincher(path: str) -> dict:
    from bioforge.evidence.readers import fincher
    return {"adata": fincher.read_dge(path)}


@register("ingest.plass", outputs=["adata"],
          description="Load the Plass 2018 atlas matrix as AnnData (dd_Smed_v6).")
def ingest_plass(path: str) -> dict:
    from bioforge.evidence.readers import plass
    return {"adata": plass.read_plass_matrix(path)}


@register("ingest.king", outputs=["tf_catalog", "rnai", "correlations", "atlas"],
          description="Load King 2024 supplementary xlsx files (mmc4/5/6/7).")
def ingest_king(catalog_path: str = "", rnai_path: str = "",
                correlations_path: str = "", atlas_path: str = "") -> dict:
    from bioforge.evidence.readers import king
    out: dict[str, Any] = {}
    if catalog_path:
        out["tf_catalog"] = king.read_tf_catalog(catalog_path)
    if rnai_path:
        out["rnai"] = king.read_rnai_table(rnai_path)
    if correlations_path:
        out["correlations"] = king.read_correlations(correlations_path)
    if atlas_path:
        out["atlas"] = king.read_king_atlas(atlas_path)
    return out


@register("evidence.load_bridge_csv", outputs=["bridge"],
          description="Load a v4<->v6 bridge CSV into a BridgeTable.")
def load_bridge_csv(path: str) -> dict:
    from bioforge.evidence import load_bridge
    return {"bridge": load_bridge(path)}


@register("evidence.score_per_atlas", outputs=["records_per_atlas"],
          description="Score TF candidates per atlas (called inside a loop: over atlases).")
def score_per_atlas(adata: "ad.AnnData",
                    tf_catalog: "pd.DataFrame | None" = None,
                    bridge: "Any | None" = None,
                    min_log2fc: float = 1.5,
                    max_pval: float = 1e-3,
                    atlas_name: str = "") -> dict:
    '''Compute expression + specificity evidence for one atlas's AnnData.

    `tf_catalog` is the King mmc4 DataFrame (gene_id + is_tf flags) used to
    restrict the search to putative TFs. `bridge` is a BridgeTable so v4
    AnnData can be projected onto the v6 catalog for scoring.

    Returns {"records": List[EvidenceRecord]}.
    '''
    import scanpy as sc
    from bioforge.evidence.schema import EvidenceRecord, EvidenceSource

    rng = np.random.default_rng(0)  # deterministic for reproducibility
    scores: list = list(adata.var_names)

    # Filter to TFs when available
    if tf_catalog is not None and "gene_id" in tf_catalog.columns:
        tf_ids = set(tf_catalog.loc[tf_catalog["is_tf"], "gene_id"])
        if bridge is not None:
            # bridge v4 -> v6 then intersect with tf_ids
            v6_of = {v: bridge.v4_to_v6(v) for v in adata.var_names}
            keep = [v for v, v6 in v6_of.items() if v6 in tf_ids]
        else:
            # already v6
            keep = [v for v in adata.var_names if v in tf_ids]
        scores = keep or list(adata.var_names)

    records: list[EvidenceRecord] = []
    # Determine cluster column (leiden is set by Layer 8A cluster step)
    cluster_key = "leiden" if "leiden" in adata.obs else None
    if cluster_key:
        sc.tl.rank_genes_groups(adata, cluster_key, method="wilcoxon")
        result = adata.uns["rank_genes_groups"]
        for gene in scores:
            mask = result["names"] == gene
            if not mask.any():
                continue
            log2fc = result["logfoldchanges"][mask][0]
            pval = result["pvals_adj"][mask][0]
            if log2fc < min_log2fc or pval > max_pval:
                continue
            r = EvidenceRecord(
                gene_id=gene,
                gene_name=(adata.var.loc[gene, "gene_name"] if "gene_name" in adata.var else gene),
            )
            r.add_score(EvidenceSource.EXPRESSION, float(min(1.0, log2fc / 5.0)),
                        note=f"log2FC={log2fc:.2f},pval={pval:.2g},atlas={atlas_name}")
            # Specificity: complement of a Shannon entropy across clusters
            try:
                expr_per_cluster = np.asarray([
                    (adata.X[adata.obs[cluster_key] == c, :][:, adata.var_names.get_loc(gene)] > 0)
                    .sum() for c in np.unique(adata.obs[cluster_key])]).astype(float)
                probs = expr_per_cluster / (expr_per_cluster.sum() + 1e-9)
                entropy = -np.sum(probs * np.log(probs + 1e-9))
                specificity = float(np.exp(-entropy))
            except Exception:
                specificity = 0.5
            r.add_score(EvidenceSource.SPECIFICITY, float(min(1.0, specificity)),
                        note=f"atlas={atlas_name}")
            records.append(r)
    # If no ranks (e.g. cluster step skipped), still return empty records list
    return {"records": records, "atlas_name": atlas_name}


@register("evidence.combine_scores", outputs=["records"],
          description="Combine per-atlas record lists into one merged list (reproducibility stream).")
def combine_scores(records_per_atlas: list, **kwargs) -> dict:
    '''records_per_atlas: list of {"records": List[EvidenceRecord], "atlas_name": str}.'''
    merged: dict[str, EvidenceRecord] = {}
    atlas_hits: dict[str, set[str]] = {}
    for entry in records_per_atlas:
        records = entry["records"] if isinstance(entry, dict) else entry
        atlas_name = entry.get("atlas_name", "?") if isinstance(entry, dict) else "?"
        for r in records:
            if r.gene_id not in merged:
                merged[r.gene_id] = EvidenceRecord(gene_id=r.gene_id, gene_name=r.gene_name)
            existing = merged[r.gene_id]
            for src, val in r.scores.items():
                val_max = max(existing.scores.get(src, 0.0), val)
                existing.scores[src] = val_max
                if src == EvidenceSource.EXPRESSION:
                    merged_notes = existing.notes.get(src, "") + f",{atlas_name}"
                    existing.notes[src] = merged_notes
            atlas_hits.setdefault(r.gene_id, set()).add(atlas_name)
    # Reproducibility score: how many atlases supported this gene
    for g, r in merged.items():
        n_atlases = len(atlas_hits[g])
        r.add_score(EvidenceSource.REPRODUCIBILITY, float(n_atlases / 3.0),
                    note=f"atlases={sorted(atlas_hits[g])}")
    return {"records": list(merged.values()), "atlas_membership": {g: sorted(v) for g, v in atlas_hits.items()}}


@register("evidence.add_rnai_stream", inputs=["records"], outputs=["records"],
          description="Score the RNAi stream by cross-referencing King mmc5.")
def add_rnai_stream(records: list[EvidenceRecord], rnai_table: "pd.DataFrame | None" = None,
                    bridge: "Any | None" = None) -> dict:
    if rnai_table is None or "fstf_rnai" not in rnai_table.columns:
        return {"records": records}
    rnai_targets = set(rnai_table["fstf_rnai"].astype(str))
    for r in records:
        name = (r.gene_name or "")
        in_table = name in rnai_targets or r.gene_id in rnai_targets
        if bridge is not None:
            v6 = bridge.v4_to_v6(r.gene_id) or r.gene_id
            in_table = in_table or v6 in rnai_targets
        r.add_score(EvidenceSource.RNai, 1.0 if in_table else 0.0,
                    note=f"in_mmc5={in_table}")
    return {"records": records}


@register("evidence.add_correlation_stream", outputs=["records"],
          description="Score the correlation stream from King mmc6.")
def add_correlation_stream(records: list[EvidenceRecord],
                           correlations: "pd.DataFrame | None" = None) -> dict:
    if correlations is None or "tf1" not in correlations.columns:
        return {"records": records}
    relevant = correlations[correlations["tf1"].isin([r.gene_name for r in records if r.gene_name])] \
        if "tf1" in correlations.columns else correlations.iloc[0:0]
    for r in records:
        name = r.gene_name or ""
        hits = relevant[(relevant["tf1"] == name) | (relevant["tf2"] == name)]
        if len(hits) == 0:
            continue
        x1 = hits["x1_corr"].mean() if "x1_corr" in hits else 0.0
        g0 = hits["g0_corr"].mean() if "g0_corr" in hits else 0.0
        # Score = how much G0 corr exceeds X1 (the "maturation-fueled" signal)
        gain = max(0.0, float(g0) - float(x1))
        r.add_score(EvidenceSource.CORRELATION, min(1.0, gain * 3.0),
                    note=f"x1={x1:.2f},g0={g0:.2f}")
    return {"records": records}


@register("evidence.add_function_stream", outputs=["records"],
          description="Score the function stream via ontology.annotate_function.")
def add_function_stream(records: list[EvidenceRecord]) -> dict:
    from bioforge.evidence.ontology import annotate_function
    for r in records:
        annotate_function(r)
    return {"records": records}


@register("evidence.assign_tiers", outputs=["tiered"],
          description="Run assign_tiers over records, returning ordered (record, tier, score) tuples.")
def assign_tiers_step(records: list[EvidenceRecord]) -> dict:
    from bioforge.evidence import EvidenceScorer, assign_tiers
    scorer = EvidenceScorer()
    tiered = assign_tiers(records, scorer=scorer)
    # Return ordered records back-into a list for downstream steps
    tiered.sort(key=lambda t: t[2], reverse=True)
    return {"tiered": tiered, "records": [t[0] for t in tiered]}


@register("storage.save_anndata", outputs=["path"],
          description="Persist an AnnData to h5ad at the given path.")
def save_anndata(adata: "ad.AnnData", path: str) -> dict:
    import os
    os.makedirs(Path(path).parent, exist_ok=True)
    adata.write_h5ad(path)
    return {"path": path}
