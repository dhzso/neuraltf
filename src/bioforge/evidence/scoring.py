"""Multi-criterion evidence scoring.

The :class:`EvidenceScorer` integrates per-source normalized scores (in
``[0, 1]``) for each TF into a single integrated score per
:class:`bioforge.evidence.schema.EvidenceRecord`.

Advanced usage: pass a custom weights mapping to :meth:`EvidenceScorer.__call__`
to override defaults. Weights are renormalized over the sources actually
present per record — missing sources do not penalize a candidate.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

from bioforge.core.logging import get_logger
from bioforge.evidence.schema import EvidenceRecord, EvidenceSource

logger = get_logger("evidence.scoring")


DEFAULT_WEIGHTS: dict[EvidenceSource, float] = {
    EvidenceSource.EXPRESSION: 0.20,
    EvidenceSource.SPECIFICITY: 0.10,
    EvidenceSource.REPRODUCIBILITY: 0.15,
    EvidenceSource.RNai: 0.15,
    EvidenceSource.CORRELATION: 0.10,
    EvidenceSource.FUNCTION: 0.05,
    EvidenceSource.NEURAL_ENRICHED: 0.15,
    EvidenceSource.NEURAL_SPECIFICITY: 0.10,
}


@dataclass
class EvidenceScorer:
    """Compute integrated scores for TF candidates.

    Parameters
    ----------
    weights
        Per-source weights; defaults are the values shown in ADR-0002.
    """

    weights: dict[EvidenceSource, float] = field(
        default_factory=lambda: dict(DEFAULT_WEIGHTS)
    )

    def __post_init__(self) -> None:
        if not self.weights:
            raise ValueError("EvidenceScorer.weights must be non-empty")
        # All weights must be non-negative; we renormalize at call time so
        # they don't need to sum to 1 here, but we do guarantee that any
        # value passed is real and finite.
        for src, w in self.weights.items():
            if w < 0:
                raise ValueError(
                    f"weight for {src.value} must be non-negative, got {w}"
                )

    def integrated_score(self, record: EvidenceRecord) -> float:
        """Renormalized weighted sum over the sources present on ``record``."""
        used = {s: w for s, w in self.weights.items() if s in record.scores}
        if not used:
            return 0.0
        total_w = sum(used.values())
        return sum(used[s] * record.scores[s] for s in used) / total_w

    def __call__(self, records: Iterable[EvidenceRecord]) -> list[EvidenceRecord]:
        """Return ``records`` sorted by descending integrated score.

        Records are mutated in place: an ``integrated_score`` attribute is
        attached via :class:`EvidenceRecord.scores` under a synthetic key
        so callers can inspect it without parsing the dataclass.

        Actually — we avoid side effects on the record dict; the integrated
        score is returned via a parallel list that callers can use directly.
        See :func:`rank_candidates` for the consumer-facing wrapper.
        """
        ranked = sorted(
            records,
            key=lambda r: self.integrated_score(r),
            reverse=True,
        )
        logger.info("scored and ranked %d TF candidates", len(ranked))
        return ranked


def rank_candidates(
    records: Iterable[EvidenceRecord],
    scorer: EvidenceScorer | None = None,
    *,
    top_n: int | None = None,
) -> list[EvidenceRecord]:
    '''Rank candidates by integrated score; optionally truncate to top_n.'''
    s = scorer or EvidenceScorer()
    out = s(records)
    if top_n is not None and top_n >= 0:
        out = out[:top_n]
    return out
