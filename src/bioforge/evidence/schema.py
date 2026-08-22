"""Dataclasses and enums describing evidence records.

A :class:`EvidenceRecord` is the unit that the 8B framework exchanges
between modules and ultimately attaches to a TF candidate. The framework
intentionally works with simple dataclasses (not pydantic models) so it
remains dependency-light and easy to test.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class EvidenceSource(str, Enum):
    """Evidence streams integrated by the framework.

    Streams prefixed ``NEURAL_`` are thesis-specific: they encode whether
    a TF is enriched in King-atlas neural G0 progenitor subclusters, and
    how specific that neural enrichment is (fewer subclusters = higher
    specificity). These streams are populated only when the pipeline
    targets neural fate. When running the general tissue-oblivious
    pipeline they remain absent and are silently skipped by the scorer.
    """

    EXPRESSION = "expression"
    SPECIFICITY = "specificity"
    REPRODUCIBILITY = "reproducibility"
    RNai = "rnai"
    CORRELATION = "correlation"
    NEURAL_ENRICHED = "neural_enriched"
    NEURAL_SPECIFICITY = "neural_specificity"


class ConfidenceTier(str, Enum):
    """Discrete confidence label attached to ranked candidates."""

    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


@dataclass
class EvidenceRecord:
    """One TF candidate plus its per-source scores and metadata.

    Attributes
    ----------
    gene_id
        Canonical gene identifier in the working genome build
        (e.g. ``dd_Smed_v6_15104_0_1``).
    gene_name
        Human-readable gene name if known (e.g. ``"soxB"``); may be ``None``.
    scores
        Mapping of :class:`EvidenceSource` to a normalized score in [0, 1].
        Missing sources are simply absent — :func:`integrated_score` skips
        them and renormalizes the weights automatically.
    notes
        Free-form evidence annotations (e.g. which subclusters the TF was
        enriched in, the marker gene lost on RNAi, the human homolog name).
    """

    gene_id: str
    gene_name: Optional[str] = None
    scores: dict[EvidenceSource, float] = field(default_factory=dict)
    notes: dict[EvidenceSource, str] = field(default_factory=dict)
    proof_status: Optional[str] = None

    def add_score(self, source: EvidenceSource, score: float, note: str = "") -> None:
        """Attach or overwrite a per-source score (clipped to [0, 1])."""
        if not 0.0 <= score <= 1.0:
            raise ValueError(
                f"evidence score for {source.value} must be in [0, 1], got {score}"
            )
        self.scores[source] = float(score)
        if note:
            self.notes[source] = note

    def supporting_streams(self) -> int:
        """Number of evidence sources with a non-zero score."""
        return sum(1 for s in self.scores.values() if s > 0.0)
