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


# Canonical evidence-stream order (single source of truth). Consumers that
# need to enumerate streams positionally (stat scripts, sensitivity, LOO,
# permutation nulls) import STREAM_ORDER instead of hard-coding a list, so
# the number of streams stays consistent across the codebase.
STREAM_ORDER: list[EvidenceSource] = [
    EvidenceSource.EXPRESSION,
    EvidenceSource.SPECIFICITY,
    EvidenceSource.REPRODUCIBILITY,
    EvidenceSource.RNai,
    EvidenceSource.CORRELATION,
    EvidenceSource.NEURAL_ENRICHED,
    EvidenceSource.NEURAL_SPECIFICITY,
    EvidenceSource.PEREZ_LINEAGE,
    EvidenceSource.PEREZ_INFLUENCE,
    EvidenceSource.FINCHER_BRAIN,
    EvidenceSource.CUI_TEMPORAL,
]

DEFAULT_WEIGHTS: dict[EvidenceSource, float] = {
    # 11 streams summing to 1.0. RNAi and correlation are the two label-bearing
    # circular streams (King mmc5/mmc6 is the ground truth) and carry 0.05 each
    # so no single stream dominates; Fincher brain and Cui 2023 add independent
    # single-cell evidence at 0.10 each. The scorer renormalizes over the
    # streams present per record, so the proportions matter more than the sum.
    EvidenceSource.EXPRESSION:          0.100,
    EvidenceSource.SPECIFICITY:         0.100,
    EvidenceSource.REPRODUCIBILITY:     0.100,
    EvidenceSource.RNai:                0.050,
    EvidenceSource.CORRELATION:         0.050,
    EvidenceSource.NEURAL_ENRICHED:     0.100,
    EvidenceSource.NEURAL_SPECIFICITY:  0.100,
    EvidenceSource.PEREZ_LINEAGE:       0.100,   # Perez 2025 TF lineage evidence
    EvidenceSource.PEREZ_INFLUENCE:     0.100,   # Perez 2025 ANANSE regulatory influence
    EvidenceSource.FINCHER_BRAIN:       0.100,   # Fincher 2018 BrainClustering neuronal evidence
    EvidenceSource.CUI_TEMPORAL:        0.100,   # Cui 2023 regeneration time-course evidence
}


@dataclass
class EvidenceScorer:
    """Compute integrated scores for TF candidates.

    Parameters
    ----------
    weights
        Per-source weights; defaults are the values recorded in ADR-0002.
    """

    weights: dict[EvidenceSource, float] = field(
        default_factory=lambda: dict(DEFAULT_WEIGHTS)
    )

    def __post_init__(self) -> None:
        if not self.weights:
            raise ValueError("EvidenceScorer.weights must be non-empty")
        # Weights are renormalized at call time, so they don't have to sum to
        # 1 here, but every value must be a real non-negative number.
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
        if total_w == 0:
            return 0.0
        return sum(used[s] * record.scores[s] for s in used) / total_w

    def __call__(self, records: Iterable[EvidenceRecord]) -> list[EvidenceRecord]:
        """Return ``records`` sorted by descending integrated score.

        The score is not stored on the record; :func:`rank_candidates` is the
        consumer-facing wrapper.
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
