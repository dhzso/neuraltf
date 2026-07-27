"""Confidence tier assignment.

Tiers are assigned based on two cheap signals available on every
:class:`bioforge.evidence.schema.EvidenceRecord`:

1. The number of evidence streams with a non-zero normalized score
   (``record.supporting_streams()``).
2. The integrated score computed by :class:`bioforge.evidence.scoring.EvidenceScorer`.

Thresholds are documented in ADR-0002 and overridable per project.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from bioforge.core.logging import get_logger
from bioforge.evidence.scoring import EvidenceScorer
from bioforge.evidence.schema import ConfidenceTier, EvidenceRecord

logger = get_logger("evidence.confidence")


@dataclass
class ConfidencePolicy:
    """Thresholds used by :func:`assign_tiers`.

    Attributes
    ----------
    min_streams_high
        Minimum supporting streams for a HIGH tier.
    min_score_high
        Minimum integrated score for a HIGH tier.
    min_streams_medium
        Minimum supporting streams for a MEDIUM tier.
    min_score_medium
        Minimum integrated score for a MEDIUM tier.
    """

    min_streams_high: int = 4
    min_score_high: float = 0.60
    min_streams_medium: int = 3
    min_score_medium: float = 0.35


def assign_tiers(
    records: Iterable[EvidenceRecord],
    scorer: EvidenceScorer | None = None,
    policy: ConfidencePolicy | None = None,
) -> list[tuple[EvidenceRecord, ConfidenceTier, float]]:
    """Return ``(record, tier, integrated_score)`` tuples.

    Tiers are assigned greedily: HIGH takes precedence over MEDIUM over LOW.
    Records with zero supporting streams are always LOW regardless of score
    (which is also zero).
    """
    s = scorer or EvidenceScorer()
    p = policy or ConfidencePolicy()
    out: list[tuple[EvidenceRecord, ConfidenceTier, float]] = []
    for r in records:
        score = s.integrated_score(r)
        streams = r.supporting_streams()
        if streams >= p.min_streams_high and score >= p.min_score_high:
            tier = ConfidenceTier.HIGH
        elif streams >= p.min_streams_medium and score >= p.min_score_medium:
            tier = ConfidenceTier.MEDIUM
        else:
            tier = ConfidenceTier.LOW
        out.append((r, tier, score))
    logger.info("assigned tiers to %d records", len(out))
    return out
