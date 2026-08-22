"""Unit tests for bioforge.projects.neuraltf.planmine & prioritize."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from bioforge.projects.neuraltf.prioritize import (
    attach_v4,
    assign_tracks,
    compute_composite,
    extract_gene_symbol,
    go_term_flags,
    map_v6_to_v4,
    merge_annotations,
    prepare_candidates,
    rnai_marker_notes,
    select_top,
    summarize_annotations,
)
from bioforge.projects.neuraltf.planmine import (
    domain_short_name_is_dna_binding,
    go_term_flags as pgo_term_flags,
)


# ---------------------------------------------------------------------------
# planmine.domain_short_name_is_dna_binding
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name,expected", [
    ("Homeobox_dom", True),
    ("bHLH_dom", True),
    ("Znf_C2H2", True),
    ("POU", True),
    ("TF_fork_head", True),
    ("TF_T-box", True),
    ("EGF-like", False),
    ("Immunoglobulin_sub", False),
    ("S-100", False),
    (None, False),
    ("", False),
    ("nan", False),
])
def test_dna_binding_classifier(name, expected):
    assert domain_short_name_is_dna_binding(name) is expected


@pytest.mark.parametrize("term,neural,tf", [
    ("neuron differentiation", True, False),
    ("brain development", True, False),
    ("transcription factor activity, sequence-specific DNA binding", False, True),
    ("DNA binding", False, True),
    ("mitochondrial translation", False, False),
])
def test_go_term_flags(term, neural, tf):
    n, t = pgo_term_flags(term)
    assert n is neural
    assert t is tf


# ---------------------------------------------------------------------------
# v6 <-> v4 mapping
# ---------------------------------------------------------------------------
def test_map_v6_to_v4_unique_and_ambiguous():
    bridge = pd.DataFrame({
        "v6_id": ["dd_Smed_v6_1_0_1", "dd_Smed_v6_1_0_1", "dd_Smed_v6_2_0_1"],
        "v4_id": ["dd_Smed_v4_5_0_1", "dd_Smed_v4_9_0_1", "dd_Smed_v4_7_0_1"],
    })
    m = map_v6_to_v4(bridge)
    sets = dict(zip(m["v6_id"], m["mapping_flag"]))
    assert sets["dd_Smed_v6_2_0_1"] == "unique"
    assert sets["dd_Smed_v6_1_0_1"] == "ambiguous"  # two v4 ids -> flagged, blank
    assert m.loc[m["v6_id"] == "dd_Smed_v6_1_0_1", "v4_id"].iloc[0] == ""


def test_attach_v4_keeps_unmapped_blank():
    rank = pd.DataFrame({"gene_id": ["dd_Smed_v6_1_0_1", "dd_Smed_v6_99_0_1"]})
    mapping = pd.DataFrame({
        "v6_id": ["dd_Smed_v6_1_0_1"], "v4_id": ["dd_Smed_v4_5_0_1"],
        "mapping_flag": ["unique"],
    })
    out = attach_v4(rank, mapping)
    assert out["gene_id_v4"].iloc[0] == "dd_Smed_v4_5_0_1"
    assert out["gene_id_v4"].iloc[1] == ""
    assert out["v4_mapping_flag"].iloc[1] == "unmapped"


# ---------------------------------------------------------------------------
# feature extraction helpers
# ---------------------------------------------------------------------------
def test_extract_gene_symbol_heuristics():
    assert extract_gene_symbol("ALX homeobox protein 1 isoform X1 [Homo sapiens]") == "ALX"
    assert extract_gene_symbol("LIM/homeobox protein Lhx6 isoform 6 [Homo sapiens]") == "Lhx6"
    assert extract_gene_symbol("") == ""


def test_summarize_annotations_pivots_long_format():
    long = pd.DataFrame({
        "gene_id_v6": ["g1"] * 4,
        "kind": ["base", "go", "domain", "domain"],
        "key": ["", "GO:0001", "PFAM", "PFAM"],
        "value": ["", "transcription factor activity, DNA binding", "Homeobox_dom", "EGF-like"],
        "namespace": ["", "molecular_function", "", ""],
        "contig_length": [100, None, None, None],
    })
    s = summarize_annotations(long).set_index("gene_id_v6")
    assert s.loc["g1", "n_domains"] == 2
    assert "Homeobox_dom" in s.loc["g1", "dna_binding_domains"]  # only DBD kept
    assert "EGF-like" not in s.loc["g1", "dna_binding_domains"]
    assert s.loc["g1", "contig_length"] == 100


def test_merge_annotations_annotation_values_win():
    ann = pd.DataFrame({
        "gene_id_v6": ["g1"], "contig_length": [200],
        "dna_binding_domains": ["Homeobox_dom"], "go_terms": ["brain development"],
        "domains_all": ["Homeobox_dom"], "go_namespaces": ["BP"],
        "go_ids": ["GO:0001"], "planmine_human_ortholog_desc": ["d"],
    })
    rank = pd.DataFrame({"gene_id": ["g1"], "integrated_score": [0.5],
                         "proof_status": ["novel_candidate"], "gene_name": ["g1"]})
    pre = prepare_candidates(rank)  # initials with default feature columns
    merged = merge_annotations(pre, ann)
    assert merged["dna_binding_domains"].iloc[0] == "Homeobox_dom"
    assert merged["contig_length"].iloc[0] == 200
    assert merged["human_ortholog"].iloc[0] == ""


# ---------------------------------------------------------------------------
# scoring + selection
# ---------------------------------------------------------------------------
def _rank(n=8):
    return pd.DataFrame({
        "gene_id": [f"g{i}" for i in range(n)],
        "gene_name": [f"g{i}" for i in range(n)],
        "integrated_score": [0.9 - 0.05 * i for i in range(n)],
        "n_streams": [5] * n,
        "proof_status": ["known_rnai_validated"] * (n // 2) +
                        ["novel_candidate"] * (n - n // 2),
    })


def test_select_top_respects_composite_ties():
    df = _rank()
    df["dna_binding_domains"] = ["Homeobox_dom"] * len(df)
    df["go_terms"] = ["brain development"] * len(df)
    df = compute_composite(df)
    a, b = assign_tracks(df)
    sel = select_top(a, 5)
    assert len(sel) == len(a)  # only 4 Track-A rows available
    assert sel["rank"].tolist() == [1, 2, 3, 4]
    # selection is stable under re-sorting: composite desc == given order
    assert sel["composite_score"].is_monotonic_decreasing


def test_composite_bonus_once_per_category():
    # many neural GO terms must not farm multiple bonuses
    df = pd.DataFrame({
        "gene_id": ["g"], "gene_name": ["g"],
        "integrated_score": [0.7], "n_streams": [5],
        "proof_status": ["novel_candidate"],
        "dna_binding_domains": ["Homeobox_dom"],
        "go_terms": ["brain development; neurogenesis; neuron differentiation; "
                     "synaptic transmission"],
        "human_ortholog": ["FOXP2 [Homo sapiens]"],
        "planmine_human_ortholog_desc": [""],
    })
    out = compute_composite(df)
    # domain + neural-GO + human-ortholog (no TF-GO terms in the list)
    expected = 0.7 + 0.05 + 0.03 + 0.02
    assert abs(out["composite_score"].iloc[0] - expected) < 1e-9
    assert out["composite_score"].iloc[0] < 1.0


def test_rnai_marker_notes_reads_mmc5_markers():
    mmc5 = pd.DataFrame({
        "fstf_rnai": ["dd123", "dd999"],
        "marker_1": ["sert", "cintillo"],
        "marker_2": [np.nan, "dd210"],
    })
    note = rnai_marker_notes(mmc5, "dd123")
    assert "sert" in note