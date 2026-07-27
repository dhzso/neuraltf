"""Cross-atlas evidence cards.

A :class:`EvidenceCard` aggregates the per-source evidence for one TF
candidate into a single user-facing object, then renders it as a markdown
fragment suitable for both `runs/<ts>/ai_summary.md` and the Streamlit UI.

The key thesis-facing concern (per ADR-0003) is that novel *untested*
candidates should be surfaced separately from RNAi-validated known TFs —
so we attach a :class:`ProofStatus` to every card. The two classes feed
into the same integrated confidence tier/model.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Iterable, Optional

from bioforge.core.logging import get_logger
from bioforge.evidence.confidence import ConfidencePolicy, assign_tiers
from bioforge.evidence.scoring import EvidenceScorer
from bioforge.evidence.schema import ConfidenceTier, EvidenceRecord, EvidenceSource

logger = get_logger("evidence.cards")


class ProofStatus(str, Enum):
    KNOWN_RNAI_VALIDATED = "known_rnai_validated"
    NOVEL_CANDIDATE = "novel_candidate"
    PRIOR_FSTF_NOT_TESTED = "prior_fstf_not_tested"


@dataclass
class EvidenceCard:
    gene_id: str
    gene_name: Optional[str]
    integrated_score: float
    tier: ConfidenceTier
    proof_status: ProofStatus
    supporting_streams: int
    per_source: dict[EvidenceSource, tuple[float, str]] = field(default_factory=dict)
    atlases_supported: set[str] = field(default_factory=set)
    suggested_followups: list[str] = field(default_factory=list)


def _classify_proof(
    record: EvidenceRecord,
    is_prior_fstf_below_threshold: bool = False,
) -> tuple[ProofStatus, list[str]]:
    rnai_score = record.scores.get(EvidenceSource.RNai, 0.0)
    if rnai_score > 0.0:
        return ProofStatus.KNOWN_RNAI_VALIDATED, [
            "Existing RNAi phenotype supports re-using for follow-up analysis",
            "Co-stain with new candidate TFs to test combinatorial codes",
        ]
    if is_prior_fstf_below_threshold:
        return ProofStatus.PRIOR_FSTF_NOT_TESTED, [
            "Prior FSTF reported in literature but no King RNAi phenotype",
            "Test by RNAi + FISH to confirm functional role",
        ]
    return ProofStatus.NOVEL_CANDIDATE, [
        "Knock down by RNAi and assay neural cell-type markers by FISH",
        "Validate with co-expression of known neural FSTFs in neoblasts",
        "Check ortholog function in vertebrate models",
    ]


def build_evidence_card(
    record: EvidenceRecord,
    scorer: EvidenceScorer | None = None,
    policy: ConfidencePolicy | None = None,
    atlases_supported: Optional[Iterable[str]] = None,
    is_prior_fstf: bool = False,
) -> EvidenceCard:
    """Construct an :class:`EvidenceCard` for a single record."""
    s = scorer or EvidenceScorer()
    score = s.integrated_score(record)
    # Determine tier via the confidence module (single-record call).
    tiered = assign_tiers([record], scorer=s, policy=policy)[0]
    tier: ConfidenceTier = tiered[1]
    proof_status, followups = _classify_proof(
        record, is_prior_fstf_below_threshold=is_prior_fstf and record.scores.get(
            EvidenceSource.RNai, 0.0) == 0.0
    )
    return EvidenceCard(
        gene_id=record.gene_id,
        gene_name=record.gene_name,
        integrated_score=score,
        tier=tier,
        proof_status=proof_status,
        supporting_streams=record.supporting_streams(),
        per_source={
            src: (val, record.notes.get(src, ""))
            for src, val in record.scores.items()
        },
        atlases_supported=set(atlases_supported or set()),
        suggested_followups=followups,
    )


def build_cards_for_records(
    records: Iterable[EvidenceRecord],
    *,
    scorer: EvidenceScorer | None = None,
    policy: ConfidencePolicy | None = None,
    atlas_membership: Optional[dict[str, set[str]]] = None,
    prior_fstf_ids: Optional[set[str]] = None,
) -> list[EvidenceCard]:
    '''Build cards for an iterable of records.

    `atlas_membership` maps gene_id → set of atlas names that supported it.
    `prior_fstf_ids` is an optional set of gene_ids that were prior FSTFs
    but had no RNAi phenotype, so the framework marks them
    PRIOR_FSTF_NOT_TESTED rather than NOVEL_CANDIDATE.
    '''
    s = scorer or EvidenceScorer()
    membership = atlas_membership or {}
    prior = prior_fstf_ids or set()
    cards: list[EvidenceCard] = []
    for r in records:
        cards.append(build_evidence_card(
            r, scorer=s, policy=policy,
            atlases_supported=membership.get(r.gene_id, set()),
            is_prior_fstf=r.gene_id in prior,
        ))
    cards.sort(key=lambda c: c.integrated_score, reverse=True)
    return cards


def render_card_markdown(card: EvidenceCard) -> str:
    """Render one card as a markdown fragment."""
    lines = [
        f"## {card.gene_name or card.gene_id}",
        "",
        f"- **Gene ID:** `{card.gene_id}`",
        f"- **Integrated score:** {card.integrated_score:.3f}",
        f"- **Confidence tier:** {card.tier.value}",
        f"- **Proof status:** {card.proof_status.value}",
        f"- **Supporting streams:** {card.supporting_streams}",
    ]
    if card.atlases_supported:
        atlases = ", ".join(sorted(card.atlases_supported))
        lines.append(f"- **Atlases supported:** {atlases}")
    if card.per_source:
        lines.append("")
        lines.append("### Per-source evidence")
        for src in EvidenceSource:
            if src in card.per_source:
                score, note = card.per_source[src]
                line = f"- {src.value}: {score:.2f}"
                if note:
                    line += f"  -  {note}"
                lines.append(line)
    if card.suggested_followups:
        lines.append("")
        lines.append("### Suggested follow-ups")
        for sugg in card.suggested_followups:
            lines.append(f"- {sugg}")
    return "\n".join(lines)


def render_cards_markdown(cards: list[EvidenceCard]) -> str:
    '''Render multiple cards as a single markdown report.'''
    return "\n\n---\n\n".join(render_card_markdown(c) for c in cards)
