"""Ontology mapping (stub).

Stubbed from ADR-0002; concrete functional-category and GO mapping will be
added when NeuralTF (Layer 9) provides per-tissue curated category lists.
The function signature is stable now so callers can plug in later.
"""
from __future__ import annotations

from typing import Optional

from bioforge.core.logging import get_logger
from bioforge.evidence.schema import EvidenceRecord, EvidenceSource

logger = get_logger("evidence.ontology")


# A minimal hard-coded mapping of well-known planarian TFs to functional
# categories. Real implementation will read a project-shipped CSV.
_KNOWN_CATEGORIES: dict[str, str] = {
    "soxB": "neural",
    "pou4l-1": "neural",
    "coe": "neural",
    "pitx": "neural",
    "lhx1/5-1": "neural",
    "phox2a": "neural",
    "pou4f3": "neural",
    "tbx2/3b": "neural",
    "myoD": "muscle",
    "nk4": "muscle",
    "nkx1-1": "muscle",
    "gata4/5/6-1": "intestine",
    "hnf4": "intestine",
    "prox-1": "intestine",
    "foxA": "parenchyma",
    "IRX1": "parenchyma",
    "GCM2": "parenchyma",
}


def annotate_function(
    record: EvidenceRecord, name_map: Optional[dict[str, str]] = None
) -> str:
    """Return the canonical functional category for `record.gene_name`.

    If the gene_name is unknown the function returns ``"unknown"`` and
    attaches a *zero* :class:`EvidenceSource.FUNCTION` score to the record
    when `name_map` indicates no homolog, or ``1.0`` when one is found.
    """
    name = record.gene_name or ""
    table = name_map if name_map is not None else _KNOWN_CATEGORIES
    category = table.get(name, "unknown")
    record.add_score(
        EvidenceSource.FUNCTION,
        score=1.0 if category != "unknown" else 0.0,
        note=f"category={category}",
    )
    return category
