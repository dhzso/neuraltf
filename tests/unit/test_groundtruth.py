"""Unit tests for the corrected King-2024 ground-truth labels.

2026-09-11 ground-truth correction: mmc5 lists ALL inhibited TFs
("All Transcription Factors Inhibited"); the distributed xlsx copy lost
the red/green phenotype font encoding, so `tested` means SCREENED.
These tests pin the semantics of
bioforge.evidence.groundtruth so no consumer can silently regress to
"tested == validated".
"""
from __future__ import annotations

import pytest

from bioforge.evidence.cards import (
    EvidenceCard,
    ProofStatus,
    _classify_proof,
    build_evidence_card,
)
from bioforge.evidence.groundtruth import (
    MMC5_GROUND_TRUTH_NOTE,
    PHENOTYPE_CONFIRMED_NAMES,
    PHENOTYPE_CONFIRMED_SHORT,
    PHENOTYPE_CONFIRMED_V6,
    is_phenotype_confirmed,
    phenotype_status,
)
from bioforge.evidence.schema import EvidenceRecord, EvidenceSource


class TestCuratedSet:
    def test_set_is_frozen_and_nonempty(self):
        assert isinstance(PHENOTYPE_CONFIRMED_V6, frozenset)
        assert len(PHENOTYPE_CONFIRMED_V6) == 20

    def test_all_entries_are_v6_ids(self):
        for g in PHENOTYPE_CONFIRMED_V6:
            assert g.startswith("dd_Smed_v6_"), g
            assert g.endswith("_0_1"), g

    def test_short_tokens_derived_consistently(self):
        assert len(PHENOTYPE_CONFIRMED_SHORT) == len(PHENOTYPE_CONFIRMED_V6)
        for tok in PHENOTYPE_CONFIRMED_SHORT:
            assert tok.startswith("dd")
            assert tok[2:].isdigit()

    def test_known_members(self):
        # Every FISH panel TF from the paper must be in the set
        for gid in (
            "dd_Smed_v6_29211_0_1",  # PHOX2A
            "dd_Smed_v6_22163_0_1",  # UNCX
            "dd_Smed_v6_14753_0_1",  # ascl-2
            "dd_Smed_v6_8096_0_1",   # fer3l-1
            "dd_Smed_v6_8104_0_1",   # SOX2 / soxB1-2
            "dd_Smed_v6_6470_0_1",   # Tbx2/3b (via build_king_atlas manual map)
            "dd_Smed_v6_11693_0_1",  # Tbx2/3c
            "dd_Smed_v6_9061_0_1",   # Post-2b
        ):
            assert gid in PHENOTYPE_CONFIRMED_V6, gid

    def test_screened_without_phenotype_excluded(self):
        # dd17143 (TBX2) was screened but shows no published phenotype
        assert "dd_Smed_v6_17143_0_1" not in PHENOTYPE_CONFIRMED_V6


class TestIsPhenotypeConfirmed:
    def test_exact_v6_match(self):
        assert is_phenotype_confirmed("dd_Smed_v6_29211_0_1")

    def test_numeric_token_match(self):
        # v4-dialect ids share the numeric token (v4 tables are isoform-blind)
        assert is_phenotype_confirmed("dd_Smed_v4_29211_0_1")

    def test_v6_sibling_isoform_does_not_inherit_flag(self):
        # 2026-09-06 audit rule: distinct v6 isoforms stay distinct. The
        # curated entry is dd_Smed_v6_10911_0_1; its numeric sibling
        # dd_Smed_v6_10911_0_2 must NOT inherit the phenotype flag.
        assert "dd_Smed_v6_10911_0_1" in PHENOTYPE_CONFIRMED_V6
        assert not is_phenotype_confirmed("dd_Smed_v6_10911_0_2")

    def test_negative(self):
        assert not is_phenotype_confirmed("dd_Smed_v6_17143_0_1")
        assert not is_phenotype_confirmed("dd_Smed_v6_38342_0_1")  # top-1 tested, screened only
        assert not is_phenotype_confirmed("dd_Smed_v6_99999_0_1")
        assert not is_phenotype_confirmed("dd_Smed_v6_10911_0_2")  # sibling isoform

    def test_garbage_safe(self):
        assert not is_phenotype_confirmed(None)
        assert not is_phenotype_confirmed("")
        assert not is_phenotype_confirmed("nan")
        assert not is_phenotype_confirmed("dd_Smed_v6_")
        assert not is_phenotype_confirmed("pax2b")  # name forms never match


class TestPhenotypeStatus:
    def test_confirmed(self):
        assert phenotype_status(1.0, "dd_Smed_v6_29211_0_1") == "phenotype_confirmed"

    def test_screened_no_phenotype(self):
        assert phenotype_status(1.0, "dd_Smed_v6_38342_0_1") == "screened_no_phenotype"
        assert phenotype_status(1.0, "dd_Smed_v6_17143_0_1") == "screened_no_phenotype"

    def test_not_screened(self):
        assert phenotype_status(0.0, "dd_Smed_v6_5882_0_1") == "not_screened"
        assert phenotype_status(None, "dd_Smed_v6_5882_0_1") == "not_screened"
        # confirmed genes that were never matched to the screening list
        # (name-form gaps like GCM2/Tbx2/3b) still report their phenotype
        assert phenotype_status(0.0, "dd_Smed_v6_6470_0_1") == "not_screened"


class TestProofStatusSemantics:
    def _record(self, gene_id, rnai):
        r = EvidenceRecord(gene_id=gene_id)
        r.add_score(EvidenceSource.RNai, rnai)
        return r

    def test_tested_means_screened_not_validated(self):
        # The enum docstring is the contract: TESTED == screened.
        doc = ProofStatus.TESTED.__doc__.lower()
        assert "screened" in doc
        assert "phenotype not implied" in doc

    def test_classify_confirmed_gets_flagged_followups(self):
        rec = self._record("dd_Smed_v6_29211_0_1", 1.0)
        status, followups = _classify_proof(rec)
        assert status is ProofStatus.TESTED
        assert any("FISH-confirmed" in f for f in followups)
        card = build_evidence_card(rec)
        assert card.phenotype_confirmed is True
        assert rec.phenotype_confirmed is True

    def test_classify_screened_only_gets_caveat_followups(self):
        rec = self._record("dd_Smed_v6_38342_0_1", 1.0)
        status, followups = _classify_proof(rec)
        assert status is ProofStatus.TESTED
        assert any("no published phenotype" in f for f in followups)
        card = build_evidence_card(rec)
        assert card.phenotype_confirmed is False

    def test_record_carries_phenotype_flag(self):
        rec = self._record("dd_Smed_v6_22163_0_1", 1.0)
        build_evidence_card(rec)
        assert rec.phenotype_confirmed is True
        rec2 = self._record("dd_Smed_v6_12722_0_1", 1.0)
        build_evidence_card(rec2)
        assert rec2.phenotype_confirmed is False

    def test_not_tested_default_flag_false(self):
        rec = self._record("dd_Smed_v6_5882_0_1", 0.0)
        card = build_evidence_card(rec)
        assert card.proof_status is ProofStatus.NOT_TESTED
        assert card.phenotype_confirmed is False

    def test_card_markdown_shows_confirmed_marker(self):
        from bioforge.evidence.cards import render_card_markdown
        rec = self._record("dd_Smed_v6_29211_0_1", 1.0)
        card = build_evidence_card(rec)
        md = render_card_markdown(card)
        assert "(phenotype-confirmed †)" in md


class TestGroundTruthNote:
    def test_note_mentions_screening_title(self):
        assert "All Transcription Factors Inhibited" in MMC5_GROUND_TRUTH_NOTE
        assert "SCREENED" in MMC5_GROUND_TRUTH_NOTE

    def test_names_aliases_exist(self):
        # name-form aliases for figure-labeled genes
        for name in ("phox2a", "uncx", "ascl-2", "fer3l-1", "tbx2/3b", "sox2"):
            assert name in PHENOTYPE_CONFIRMED_NAMES


class TestEvidenceRecordField:
    def test_phenotype_confirmed_defaults_false(self):
        r = EvidenceRecord(gene_id="dd_Smed_v6_1_0_1")
        assert r.phenotype_confirmed is False


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
