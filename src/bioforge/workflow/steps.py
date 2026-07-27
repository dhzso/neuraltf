"""Built-in workflow steps wiring Layers 8A / 8B to the Layer 7 registry.

This module is imported lazily by the CLI / UI so the framework's omics
dependencies don't have to be importable when only the registry is used.
Every callable here has a small docstring so an LLM-driven assistant can
introspect workflows through the StepRegistry.
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from bioforge.core.logging import get_logger
from bioforge.evidence.schema import EvidenceRecord, EvidenceSource
from bioforge.ingest import ingest_dataset
from bioforge.omics.batch import run_harmony
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
