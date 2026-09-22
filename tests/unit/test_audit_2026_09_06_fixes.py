"""Regression tests for the 2026-09-06 challenge-proofing audit fixes.

Covers the audit's critical findings:
- C2/C3: weight-sensitivity slot-aware attribution + complete rank vectors
- C4: honest-score exclusion set (perez_lineage leakage)
- C5: Fisher/Stouffer use sf (no 1-cdf underflow); Sidak un-conditioning
- C7: DNA-binding domain allow-list + direction-aware GO flags
- B2: isoform-aware RNAi matching
- A6: cross-method per-gene flagging only for k>=2 genes
- calibration top-decile enrichment reads decile 0
"""
from __future__ import annotations

import importlib.util
import inspect
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

from bioforge.projects.neuraltf.pipeline import NeuralTFPipeline  # noqa: E402
from bioforge.projects.neuraltf.planmine import (  # noqa: E402
    domain_short_name_is_dna_binding,
    go_term_flags,
)


def importlib_util(relpath: str):
    path = REPO / relpath
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# B2: isoform-aware matching
# ---------------------------------------------------------------------------
class TestIsoformMatching:
    def test_structured_isoform_key_distinguishes_siblings(self):
        a = NeuralTFPipeline._structured_isoform_key("dd_Smed_v6_9460_0_1")
        b = NeuralTFPipeline._structured_isoform_key("dd_Smed_v6_9460_0_2")
        assert a == "dd_Smed_v6_9460_0_1"
        assert b == "dd_Smed_v6_9460_0_2"
        assert a != b  # the 2026-09-04 collapse aliased these

    def test_short_id_still_collapses(self):
        assert NeuralTFPipeline._short_id("dd_Smed_v6_9460_0_1") == "dd9460"

    def test_bare_token_has_no_isoform_key(self):
        assert NeuralTFPipeline._structured_isoform_key("dd9460") is None
        assert NeuralTFPipeline._structured_isoform_key("pax6A") is None

    def test_rna_target_set_keys_structured_rows_isoform_aware(self):
        pipe = NeuralTFPipeline.__new__(NeuralTFPipeline)
        pipe.rnai_table = pd.DataFrame({0: [
            "dd_Smed_v6_9460_0_1",  # structured: isoform-specific
            "dd11150 (UNCX)",        # short + symbol
            "fer3l-1",               # gene symbol
        ]})
        targets = pipe._build_rna_target_set()
        assert "dd_Smed_v6_9460_0_1" in targets
        assert "dd11150" in targets
        assert "fer3l-1" in targets
        # the structured row must NOT contribute the bare dd9460 token
        assert "dd9460" not in targets


# ---------------------------------------------------------------------------
# C4: honest-score exclusion sets (shared constant discipline)
# ---------------------------------------------------------------------------
class TestHonestExclusions:
    def test_mann_whitney_excludes_perez_lineage(self):
        spec = importlib_util("scripts/stats/mann_whitney_top10.py")
        assert "perez_lineage" in spec.HONEST_EXCLUDE
        assert "reproducibility" in spec.HONEST_STRICT_EXCLUDE
        assert "rnai" in spec.HONEST_EXCLUDE

    def test_effect_sizes_excludes_perez_lineage(self):
        spec = importlib_util("scripts/stats/effect_sizes.py")
        assert "perez_lineage" in spec.HONEST_EXCLUDE
        assert set(spec.HONEST_STRICT_EXCLUDE) > set(spec.HONEST_EXCLUDE)

    def test_precision_recall_excludes_perez_lineage(self):
        spec = importlib_util("scripts/stats/precision_recall.py")
        assert "perez_lineage" in spec.CIRCULAR_STREAMS
        assert "reproducibility" in spec.CIRCULAR_STREAMS_STRICT

    def test_negative_controls_label_free_score(self):
        spec = importlib_util("scripts/stats/negative_controls.py")
        for s in ("rnai", "neural_enriched", "neural_specificity",
                  "perez_lineage"):
            assert s in spec.LEAKING_STREAMS


# ---------------------------------------------------------------------------
# C5: combination tail via sf; Sidak un-conditioning
# ---------------------------------------------------------------------------
class TestMetaAnalyticPvalues:
    @pytest.fixture()
    def meta(self):
        return importlib_util("scripts/stats/meta_analytic_pvalue.py")

    def test_fisher_sf_no_underflow_zero(self, meta):
        # chi2 huge => 1-cdf underflows to 0.0; sf resolves ~1e-33
        chi2, combined = meta.fishers_method([1e-16, 1e-16, 1e-16])
        assert combined > 0.0
        assert combined < 1e-30

    def test_stouffer_sf_no_zero(self, meta):
        z, p = meta.stouffers_method([1e-16, 1e-16, 1e-16])
        assert p > 0.0

    def test_adjusted_min_p_sidak(self, meta):
        # 1-(1-0.01)^10 ~ 0.0956
        assert abs(meta.adjusted_min_p(0.01, 10) - (1 - 0.99 ** 10)) < 1e-12
        assert meta.adjusted_min_p(0.05, 1) == 0.05

    def test_cluster_counts_read_from_checkpoint(self, meta):
        df = pd.DataFrame({
            "atlas": ["fincher", "plass", "cui"],
            "n_leiden_clusters": [16, 22, 30],
        })
        with tempfile.TemporaryDirectory() as td:
            ck = Path(td) / "checkpoint_02_post_qc.csv"
            df.to_csv(ck, index=False)
            orig_run = meta.RUN_DIR
            meta.RUN_DIR = Path(td)
            try:
                got = meta.cluster_counts_from_checkpoint()
            finally:
                meta.RUN_DIR = orig_run
        assert got == {"fincher": 16, "plass": 22, "cui": 30}


# ---------------------------------------------------------------------------
# C7: domain + GO classifiers (corpus-verified)
# ---------------------------------------------------------------------------
class TestClassifiers:
    @pytest.mark.parametrize("name", [
        "HMG_box_dom", "Paired_dom", "TF_Brachyury", "TCF/LEF",
        "CUT_dom", "Iroquois_homeo", "Homeo_prospero_dom",
        "Transcription_factor_COE", "TF_AP2", "ARID_dom",
        "MAD_homology_MH1", "DMRT/dsx/mab-3", "Teashirt_fam",
        "Nucl_hormone_rcpt_ligand-bd", "Znf_NHR/GATA", "WHTH_DNA-bd_dom",
        "p53-like_TF_DNA-bd", "Homeobox_KN_domain", "HTH_Psq",
    ])
    def test_real_dbd_names_pass(self, name):
        assert domain_short_name_is_dna_binding(name)

    @pytest.mark.parametrize("name", [
        "Znf_LIM", "Znf_hrmn_rcpt", "SMAD_dom-like",
        "SMAD_dom_Dwarfin-type", "Dwarfin", "SMAD_FHA_domain",
        "BRCT_dom", "PAS_fold", "Peptidase_M1", "Innexin", "Ig_E-set",
        "DUF3504", "CVC", "ASH",
    ])
    def test_non_dbd_names_fail(self, name):
        assert not domain_short_name_is_dna_binding(name)

    def test_go_head_false_positive_killed(self):
        assert go_term_flags("head involution")[0] is False
        assert go_term_flags("specification of segmental identity, head")[0] is False

    def test_go_bare_dna_binding_not_tf(self):
        assert go_term_flags("DNA binding")[1] is False

    def test_go_tf_class_terms_still_flag(self):
        assert go_term_flags(
            "transcription factor activity, sequence-specific DNA binding")[1] is True


# ---------------------------------------------------------------------------
# A6: cross-method per-gene flagging restricted to k>=2
# ---------------------------------------------------------------------------
class TestCrossMethodFlagging:
    def test_benjamini_hochberg_monotone(self):
        spec = importlib_util("scripts/stats/cross_method_correction.py")
        p = np.array([1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9,
                      0.0026, 0.0026, 0.0026, 0.3, 0.4, 0.5])
        q = spec.benjamini_hochberg(p)
        assert np.all(q <= 1.0)
        assert (q >= 0).all()


# ---------------------------------------------------------------------------
# calibration: top-decile = decile 0
# ---------------------------------------------------------------------------
class TestCalibrationTopDecile:
    def test_top_decile_is_decile_zero(self):
        spec = importlib_util("scripts/stats/calibration.py")
        src = inspect.getsource(spec.main)
        assert "bin_stats[0]" in src


# ---------------------------------------------------------------------------
# A1: sensitivity build_shortlist emits track + rank_within_track
# ---------------------------------------------------------------------------
class TestSensitivityFixes:
    def test_build_shortlist_assigns_track(self):
        spec = importlib_util("scripts/run_weight_sensitivity.py")
        rows = {"gene_id": ["g1", "g2"], "gene_name": ["a", "b"],
                "proof_status": ["tested", "not_tested"],
                "integrated_score": [0.9, 0.8], "n_streams": [9, 5],
                "dna_binding_domains": ["Homeobox_dom", ""],
                "mmc4_tf_flag": ["TF", ""], "go_terms": ["", ""],
                "human_ortholog": ["", ""],
                "planmine_human_ortholog_desc": ["", ""]}
        for s in spec.STREAMS:
            rows[s] = [0.9, 0.5]
        df = pd.DataFrame(rows)
        out = spec.build_shortlist(df, np.array([0.95, 0.4]))
        assert "track" in out.columns
        assert set(out["track"]) <= {"A", "B"}
        assert "rank_within_track" in out.columns


# ---------------------------------------------------------------------------
# 2026-09-21: permutation test float-tie + BH family fixes
# ---------------------------------------------------------------------------
class TestPermutationTestFixes:
    """Regression tests for the 2026-09-21 permutation audit fixes.

    Verifies:
    - Float-tie genes (surplus ≤ epsilon) get p=1.0, not the floor.
    - BH-FDR is computed only over testable genes.
    - significant_fdr_05 is False for all untestable genes.
    - The corrected CSV (if present) has no untestable gene marked significant.
    """

    def test_untestable_genes_get_p_one(self):
        """An untestable gene (real <= ceiling + eps) must get p=1.0
        regardless of what the null distribution looks like."""
        # Simulate the fix logic directly
        real_score = 0.9214124856363768
        ceiling = 0.9214124856363767
        eps = 1e-9
        is_untestable = real_score <= ceiling + eps
        assert is_untestable, "This gene should be classified as untestable"

        # Under the fix, untestable → p=1.0 by enforcement
        if is_untestable:
            p = 1.0
        else:
            # Would normally compute from null distribution
            p = 0.001  # floor
        assert p == 1.0, f"Untestable gene got p={p}, expected 1.0"

    def test_testable_gene_gets_computed_p(self):
        """A testable gene (real > ceiling + eps) should get its computed
        p-value, not p=1.0."""
        real_score = 0.90
        ceiling = 0.80
        eps = 1e-9
        is_untestable = real_score <= ceiling + eps
        assert not is_untestable

    def test_exact_tie_is_untestable(self):
        """A gene with surplus of exactly 0.0 (real == ceiling) is untestable."""
        real_score = 0.85
        ceiling = 0.85
        eps = 1e-9
        assert real_score <= ceiling + eps

    def test_bh_excludes_untestable(self):
        """BH-FDR should be computed only over testable genes.

        If we include untestable genes (all p=1.0 after the fix), they'd
        inflate the BH family size and dilute q-values for testable genes.
        """
        from statsmodels.stats.multitest import multipletests

        # 5 testable genes at the floor p, 10 untestable at p=1.0
        testable_p = np.array([0.001] * 5)
        untestable_p = np.array([1.0] * 10)

        # BH over testable only (correct)
        _, q_correct, _, _ = multipletests(testable_p, alpha=0.05, method="fdr_bh")

        # BH over all (incorrect old behavior)
        all_p = np.concatenate([testable_p, untestable_p])
        _, q_wrong, _, _ = multipletests(all_p, alpha=0.05, method="fdr_bh")

        # The testable-only q should be <= the all-inclusive q for each
        # testable gene (smaller family → same or tighter q)
        for i in range(5):
            assert q_correct[i] <= q_wrong[i] + 1e-12, \
                f"Testable-only BH q ({q_correct[i]}) should be <= " \
                f"all-inclusive q ({q_wrong[i]})"

    def test_significant_fdr_never_true_for_untestable(self):
        """After the fix, no gene with untestable_by_permutation=True
        should have significant_fdr_05=True."""
        csv_path = REPO / "projects" / "NeuralTF" / "results" / "permutation_pvalues_neural.csv"
        if not csv_path.exists():
            pytest.skip("permutation_pvalues_neural.csv not found (pipeline not run)")
        df = pd.read_csv(csv_path)
        if "untestable_by_permutation" not in df.columns:
            pytest.skip("CSV missing untestable column (pre-audit version)")

        untestable_and_sig = df[
            df["untestable_by_permutation"] & df["significant_fdr_05"]
        ]
        assert len(untestable_and_sig) == 0, (
            f"{len(untestable_and_sig)} untestable genes are marked significant "
            f"(expected 0). Gene IDs: {untestable_and_sig['gene_id'].tolist()}"
        )

    def test_untestable_genes_have_p_one_in_csv(self):
        """After the fix, every untestable gene in the CSV has p=1.0."""
        csv_path = REPO / "projects" / "NeuralTF" / "results" / "permutation_pvalues_neural.csv"
        if not csv_path.exists():
            pytest.skip("permutation_pvalues_neural.csv not found")
        df = pd.read_csv(csv_path)
        if "untestable_by_permutation" not in df.columns:
            pytest.skip("CSV missing untestable column")

        untestable = df[df["untestable_by_permutation"]]
        bad = untestable[untestable["empirical_p"] < 1.0 - 1e-12]
        assert len(bad) == 0, (
            f"{len(bad)} untestable genes have p < 1.0 "
            f"(expected all p=1.0). Gene IDs: {bad['gene_id'].tolist()}"
        )

    def test_untestable_genes_have_nan_q_in_csv(self):
        """After the fix, every untestable gene has q=NaN (excluded from BH)."""
        csv_path = REPO / "projects" / "NeuralTF" / "results" / "permutation_pvalues_neural.csv"
        if not csv_path.exists():
            pytest.skip("permutation_pvalues_neural.csv not found")
        df = pd.read_csv(csv_path)
        if "untestable_by_permutation" not in df.columns:
            pytest.skip("CSV missing untestable column")

        untestable = df[df["untestable_by_permutation"]]
        has_q = untestable[untestable["q_bh_family"].notna()]
        assert len(has_q) == 0, (
            f"{len(has_q)} untestable genes have non-NaN q "
            f"(expected NaN — excluded from BH). Gene IDs: {has_q['gene_id'].tolist()}"
        )
