"""Confidence tier assignment.

Tiers are assigned based on cheap signals available on every
:class:`bioforge.evidence.schema.EvidenceRecord`:

1. RNAi-validated (score > 0 in :attr:`EvidenceSource.RNAi`) is always
   ``HIGH`` — an existing phenotype is the strongest possible evidence.
2. The integrated score computed by :class:`bioforge.evidence.scoring.EvidenceScorer`.
3. The number of evidence streams with a non-zero normalized score
   (``record.supporting_streams()``).

Thresholds are documented in ADR-0002 and overridable per project.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from bioforge.core.logging import get_logger
from bioforge.evidence.scoring import EvidenceScorer
from bioforge.evidence.schema import ConfidenceTier, EvidenceRecord, EvidenceSource

logger = get_logger("evidence.confidence")


@dataclass
class ConfidencePolicy:
    """Thresholds used by :func:`assign_tiers`.

    Attributes
    ----------
    min_streams_high
        Minimum supporting streams for a HIGH tier (not applied when
        RNAi is present).
    min_score_high
        Minimum integrated score for a HIGH tier.
    min_streams_medium
        Minimum supporting streams for a MEDIUM tier.
    min_score_medium
        Minimum integrated score for a MEDIUM tier.
    """

    min_streams_high: int = 3
    min_score_high: float = 0.45
    min_streams_medium: int = 2
    min_score_medium: float = 0.25


def assign_tiers(
    records: Iterable[EvidenceRecord],
    scorer: EvidenceScorer | None = None,
    policy: ConfidencePolicy | None = None,
) -> list[tuple[EvidenceRecord, ConfidenceTier, float]]:
    """Return ``(record, tier, integrated_score)`` tuples.

    Tiers are assigned greedily: RNAi-validated always HIGH, then HIGH
    thresholds, then MEDIUM, then LOW.  Records with zero supporting
    streams are always LOW regardless of score (which is also zero).
    """
    s = scorer or EvidenceScorer()
    p = policy or ConfidencePolicy()
    out: list[tuple[EvidenceRecord, ConfidenceTier, float]] = []
    for r in records:
        score = s.integrated_score(r)
        streams = r.supporting_streams()
        rnai_score = r.scores.get(EvidenceSource("rnai"), 0.0)
        if rnai_score > 0.0:
            tier = ConfidenceTier.HIGH
        elif streams >= p.min_streams_high and score >= p.min_score_high:
            tier = ConfidenceTier.HIGH
        elif streams >= p.min_streams_medium and score >= p.min_score_medium:
            tier = ConfidenceTier.MEDIUM
        else:
            tier = ConfidenceTier.LOW
        out.append((r, tier, score))
    logger.info("assigned tiers to %d records", len(out))
    return out
