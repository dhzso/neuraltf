"""Tests for bioforge.evidence.cards — cross-atlas evidence cards."""
from __future__ import annotations

import pytest

from bioforge.evidence import (
    ConfidenceTier,
    EvidenceCard,
    EvidenceRecord,
    EvidenceSource,
    ProofStatus,
    build_cards_for_records,
    build_evidence_card,
    render_card_markdown,
)


def _record_with_streams(streams):
    r = EvidenceRecord(gene_id="dd_Smed_v6_42", gene_name="someTf")
    for src, val in streams:
        r.add_score(src, val)
    return r


def test_card_classifies_known_rnai_validated() -> None:
    r = _record_with_streams([
        (EvidenceSource.EXPRESSION, 1.0),
        (EvidenceSource.SPECIFICITY, 0.8),
        (EvidenceSource.REPRODUCIBILITY, 1.0),
        (EvidenceSource.RNai, 1.0),
        (EvidenceSource.CORRELATION, 0.6),
    ])
    card = build_evidence_card(r)
    assert card.proof_status == ProofStatus.KNOWN_RNAI_VALIDATED
    assert card.integrated_score > 0.6
    assert card.tier == ConfidenceTier.HIGH
    assert "RNAi phenotype" in card.suggested_followups[0]


def test_card_classifies_novel_candidate_when_no_rnai_no_prior() -> None:
    r = _record_with_streams([
        (EvidenceSource.EXPRESSION, 0.95),
        (EvidenceSource.SPECIFICITY, 0.9),
        (EvidenceSource.REPRODUCIBILITY, 0.66),
    ])
    card = build_evidence_card(r)
    assert card.proof_status == ProofStatus.NOVEL_CANDIDATE
    # suggest wet-lab followups like RNAi / FISH
    assert any("RNAi" in s for s in card.suggested_followups)


def test_card_classifies_prior_fstf_not_tested() -> None:
    r = _record_with_streams([
        (EvidenceSource.EXPRESSION, 1.0),
        (EvidenceSource.SPECIFICITY, 0.7),
        (EvidenceSource.REPRODUCIBILITY, 0.66),
    ])
    # pass is_prior_fstf=True — no RNAi score (above), so should land in prior-not-tested
    card = build_evidence_card(r, is_prior_fstf=True)
    assert card.proof_status == ProofStatus.PRIOR_FSTF_NOT_TESTED
    assert any("RNAi" in s for s in card.suggested_followups)


def test_card_includes_atlas_membership() -> None:
    r = _record_with_streams([(EvidenceSource.EXPRESSION, 1.0),
                              (EvidenceSource.REPRODUCIBILITY, 1.0)])
    card = build_evidence_card(r, atlases_supported={"fincher", "plass", "king"})
    assert card.atlases_supported == {"fincher", "plass", "king"}


def test_build_cards_for_records_sorts_by_score_desc() -> None:
    r_lo = _record_with_streams([(EvidenceSource.EXPRESSION, 0.2)])
    r_hi = _record_with_streams([(EvidenceSource.EXPRESSION, 1.0)])
    cards = build_cards_for_records(
        [r_lo, r_hi],
        atlas_membership={r_hi.gene_id: {"king"}, r_lo.gene_id: {"fincher"}},
    )
    assert cards[0].integrated_score > cards[1].integrated_score


def test_render_card_markdown_contains_key_sections() -> None:
    r = _record_with_streams([
        (EvidenceSource.EXPRESSION, 0.8),
        (EvidenceSource.SPECIFICITY, 0.5),
    ])
    card = build_evidence_card(r, atlases_supported={"king", "plass"})
    md = render_card_markdown(card)
    assert "## someTf" in md
    assert "Per-source evidence" in md
    assert "Atlases supported" in md
    assert "Suggested follow-ups" in md
    assert "expression:" in md


def test_render_cards_markdown_joins_multiple() -> None:
    from bioforge.evidence.cards import render_cards_markdown
    r1 = _record_with_streams([(EvidenceSource.EXPRESSION, 0.5)])
    r2 = _record_with_streams([(EvidenceSource.EXPRESSION, 0.8)])
    cards = build_cards_for_records([r1, r2])
    md = render_cards_markdown(cards)
    assert md.count("---") == 1  # one separator between two cards
