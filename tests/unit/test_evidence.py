"""Unit tests for bioforge.evidence (Layer 8B Evidence Integration Framework)."""
from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from bioforge.evidence import (
    AtlasHarmonizer,
    BridgeTable,
    ConfidencePolicy,
    ConfidenceTier,
    DEFAULT_WEIGHTS,
    EvidenceRecord,
    EvidenceScorer,
    EvidenceSource,
    assign_tiers,
    build_bridge_from_names,
    load_bridge,
    rank_candidates,
)
from bioforge.evidence.ontology import annotate_function


# ---------------------------------------------------------------------------
# schema.EvidenceRecord
# ---------------------------------------------------------------------------


def test_evidence_record_add_score_clips_to_unit_interval() -> None:
    r = EvidenceRecord(gene_id="dd_Smed_v6_1", gene_name="foo")
    r.add_score(EvidenceSource.EXPRESSION, 0.5)
    assert r.scores[EvidenceSource.EXPRESSION] == 0.5
    r.add_score(EvidenceSource.EXPRESSION, 1.0)
    assert r.scores[EvidenceSource.EXPRESSION] == 1.0


def test_evidence_record_add_score_rejects_out_of_range() -> None:
    r = EvidenceRecord(gene_id="x")
    with pytest.raises(ValueError):
        r.add_score(EvidenceSource.RNai, 1.5)
    with pytest.raises(ValueError):
        r.add_score(EvidenceSource.RNai, -0.1)


def test_evidence_record_supporting_streams_counts_only_nonzero() -> None:
    r = EvidenceRecord(gene_id="x")
    r.add_score(EvidenceSource.EXPRESSION, 0.0)
    r.add_score(EvidenceSource.SPECIFICITY, 0.7)
    r.add_score(EvidenceSource.RNai, 1.0)
    assert r.supporting_streams() == 2


# ---------------------------------------------------------------------------
# gene_mapping.BridgeTable / load_bridge / build_bridge_from_names
# ---------------------------------------------------------------------------


def _bridge_csv_text() -> str:
    return (
        "gene_name,v6_id,v4_id\n"
        "soxB,dd_Smed_v6_15104_0_1,dd_Smed_v4_6001_0_1\n"
        "myoD,dd_Smed_v6_22001_0_1,\n"
        "hnf4,,dd_Smed_v4_9001_0_1\n"
    )


def test_bridge_table_load_round_trip(tmp_path: Path) -> None:
    p = tmp_path / "bridge.csv"
    p.write_text(_bridge_csv_text())
    bt = load_bridge(p)
    assert bt.n_rows == 3
    assert bt.n_bridged == 1  # only soxB is present in both
    assert bt.v6_to_v4("dd_Smed_v6_15104_0_1") == "dd_Smed_v4_6001_0_1"
    assert bt.v6_to_v4("dd_Smed_v6_22001_0_1") is None
    assert bt.v4_to_v6("dd_Smed_v4_9001_0_1") is None
    assert bt.v6_to_name("dd_Smed_v6_15104_0_1") == "soxB"
    assert bt.v4_to_name("dd_Smed_v4_6001_0_1") == "soxB"


def test_bridge_table_missing_columns_raises() -> None:
    with pytest.raises(ValueError):
        BridgeTable(df=pd.DataFrame({"gene_name": ["a"], "v6_id": ["x"]}))


def test_build_bridge_from_names_outer_join() -> None:
    v6 = pd.DataFrame({
        "gene_name": ["soxB", "myoD", "unique_in_v6"],
        "v6_id": ["v6_1", "v6_2", "v6_3"],
    })
    v4 = pd.DataFrame({
        "gene_name": ["soxB", "myoD", "unique_in_v4"],
        "v4_id": ["v4_1", "v4_2", "v4_3"],
    })
    bt = build_bridge_from_names(v6, v4)
    assert bt.n_rows == 4
    assert bt.n_bridged == 2
    assert "unique_in_v6" in set(bt.df["gene_name"])
    assert "unique_in_v4" in set(bt.df["gene_name"])


def test_build_bridge_from_names_drops_duplicate_names() -> None:
    v6 = pd.DataFrame({
        "gene_name": ["soxB", "soxB", "myoD"],
        "v6_id": ["v6_a", "v6_b", "v6_c"],
    })
    v4 = pd.DataFrame({
        "gene_name": ["myoD"],
        "v4_id": ["v4_c"],
    })
    bt = build_bridge_from_names(v6, v4)
    # soxB dropped (duplicate in v6), myoD retained in both
    assert "myoD" in set(bt.df["gene_name"])
    assert "soxB" not in set(bt.df["gene_name"])
    assert bt.n_bridged == 1


# ---------------------------------------------------------------------------
# harmonization.AtlasHarmonizer
# ---------------------------------------------------------------------------


def test_atlas_harmonizer_with_default_atlases_known() -> None:
    h = AtlasHarmonizer.with_default_mappings()
    h.add_atlas("fincher", {"c1": "neuron", "c2": "muscle"})
    assert h.map("fincher", "c1") == "neuron"
    assert h.map("fincher", "unknown_label") == "unknown"


def test_atlas_harmonizer_unknown_atlas_raises() -> None:
    h = AtlasHarmonizer.with_default_mappings()
    with pytest.raises(KeyError):
        h.map("not_an_atlas", "whatever")


# ---------------------------------------------------------------------------
# scoring.EvidenceScorer / rank_candidates
# ---------------------------------------------------------------------------


def test_default_weights_sum_near_one() -> None:
    assert abs(sum(DEFAULT_WEIGHTS.values()) - 1.0) < 1e-9


def test_evidence_scorer_renormalizes_missing_sources() -> None:
    r = EvidenceRecord(gene_id="x")
    r.add_score(EvidenceSource.EXPRESSION, 1.0)  # weight 0.20
    r.add_score(EvidenceSource.RNai, 0.0)        # weight 0.15
    s = EvidenceScorer()
    # Renormalized weights: 0.20+0.15 = 0.35 -> expression gets 0.20/0.35
    # Integrated: 1.0 * (0.20/0.35) + 0.0 = 0.5714...
    assert abs(s.integrated_score(r) - (0.20 / 0.35)) < 1e-9


def test_evidence_scorer_no_scores_returns_zero() -> None:
    s = EvidenceScorer()
    assert s.integrated_score(EvidenceRecord(gene_id="x")) == 0.0


def test_evidence_scorer_rejects_negative_weight() -> None:
    with pytest.raises(ValueError):
        EvidenceScorer(weights={EvidenceSource.RNai: -0.1})


def test_rank_candidates_returns_sorted_descending() -> None:
    r_lo = EvidenceRecord(gene_id="lo")
    r_lo.add_score(EvidenceSource.EXPRESSION, 0.2)
    r_hi = EvidenceRecord(gene_id="hi")
    r_hi.add_score(EvidenceSource.EXPRESSION, 1.0)
    ranked = rank_candidates([r_lo, r_hi])
    assert [r.gene_id for r in ranked] == ["hi", "lo"]


def test_rank_candidates_top_n_truncates() -> None:
    rs = []
    for i in range(5):
        r = EvidenceRecord(gene_id=f"g{i}")
        r.add_score(EvidenceSource.EXPRESSION, i / 5)
        rs.append(r)
    top = rank_candidates(rs, top_n=2)
    assert [r.gene_id for r in top] == ["g4", "g3"]


# ---------------------------------------------------------------------------
# confidence.assign_tiers / ConfidencePolicy
# ---------------------------------------------------------------------------


def _record_with_streams(scores: list[tuple[EvidenceSource, float]]) -> EvidenceRecord:
    r = EvidenceRecord(gene_id="x")
    for src, val in scores:
        r.add_score(src, val)
    return r


def test_assign_tiers_high_band() -> None:
    streams = [
        (EvidenceSource.EXPRESSION, 1.0),
        (EvidenceSource.SPECIFICITY, 1.0),
        (EvidenceSource.REPRODUCIBILITY, 1.0),
        (EvidenceSource.RNai, 1.0),
    ]
    r = _record_with_streams(streams)
    res = assign_tiers([r])
    assert res[0][1] == ConfidenceTier.HIGH
    assert res[0][2] > 0.6


def test_assign_tiers_medium_band() -> None:
    streams = [
        (EvidenceSource.EXPRESSION, 0.5),
        (EvidenceSource.SPECIFICITY, 0.5),
    ]
    r = _record_with_streams(streams)
    res = assign_tiers([r])
    assert res[0][1] == ConfidenceTier.MEDIUM


def test_assign_tiers_low_band_when_streams_below_policy() -> None:
    streams = [(EvidenceSource.EXPRESSION, 0.2)]
    r = _record_with_streams(streams)
    res = assign_tiers([r])
    assert res[0][1] == ConfidenceTier.LOW


def test_assign_tiers_zero_streams_is_low() -> None:
    r = EvidenceRecord(gene_id="x")
    res = assign_tiers([r])
    assert res[0][1] == ConfidenceTier.LOW
    assert res[0][2] == 0.0


def test_confidence_policy_overrides_thresholds() -> None:
    streams = [(EvidenceSource.EXPRESSION, 1.0)]
    r = _record_with_streams(streams)
    strict = ConfidencePolicy(min_streams_high=1, min_score_high=0.5)
    res = assign_tiers([r], policy=strict)
    # With one source_integrated_score will be 1.0 -> HIGH under strict policy
    assert res[0][1] == ConfidenceTier.HIGH


# ---------------------------------------------------------------------------
# ontology.annotate_function
# ---------------------------------------------------------------------------


def test_annotate_function_known_tf_attaches_unit_score() -> None:
    r = EvidenceRecord(gene_id="dd_Smed_v6_99", gene_name="soxB")
    cat = annotate_function(r)
    assert cat == "neural"
    assert r.scores[EvidenceSource.FUNCTION] == 1.0


def test_annotate_function_unknown_tf_attaches_zero_score() -> None:
    r = EvidenceRecord(gene_id="dd_Smed_v6_42", gene_name="brandNewTF")
    cat = annotate_function(r)
    assert cat == "unknown"
    assert r.scores[EvidenceSource.FUNCTION] == 0.0


def test_annotate_function_custom_name_map_overrides_defaults() -> None:
    r = EvidenceRecord(gene_id="x", gene_name="customTF")
    cat = annotate_function(r, name_map={"customTF": "neuronal_progenitor"})
    assert cat == "neuronal_progenitor"
    assert r.scores[EvidenceSource.FUNCTION] == 1.0
