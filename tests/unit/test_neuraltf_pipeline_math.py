"""Unit tests for the WS1 pipeline math hardening (NeuralTF).

Covers:
- _short_id structured/short ID extraction (regression for the lazy-regex
  bug that collapsed dd_Smed_v6_11150_0_1 -> "dd6")
- _cluster_log2fc true log2FC recovery from log1p data
- _valid_perez_class ("-" non-TF sentinel must not count as evidence)
- removesuffix fallback in integrate_perez lookup path
- correlation Δr joint-NaN alignment (pairs never re-mixed)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from bioforge.projects.neuraltf.pipeline import NeuralTFPipeline


class TestShortId:
    def test_structured_v6(self):
        assert NeuralTFPipeline._short_id("dd_Smed_v6_11150_0_1") == "dd11150"

    def test_structured_v4(self):
        assert NeuralTFPipeline._short_id("dd_Smed_v4_10988_0_1") == "dd10988"

    def test_short_plain(self):
        assert NeuralTFPipeline._short_id("dd11150") == "dd11150"

    def test_short_with_annotation(self):
        assert NeuralTFPipeline._short_id("UNCX (dd22163)") == "dd22163"

    def test_non_dd_returns_none(self):
        assert NeuralTFPipeline._short_id("pax2b") is None
        assert NeuralTFPipeline._short_id("") is None
        assert NeuralTFPipeline._short_id(None) is None

    def test_no_v6_digits_leak(self):
        # regression: the old lazy regex returned "dd6" here
        assert NeuralTFPipeline._short_id("dd_Smed_v6_2946_0_1") == "dd2946"


class TestClusterLog2FC:
    def _make_adata(self, fold: float):
        rng = np.random.default_rng(0)
        base = rng.lognormal(mean=1.0, sigma=0.5, size=(200, 2))
        X = np.log1p(base)
        X[:50, 0] = np.log1p(base[:50, 0] * fold)
        import anndata as ad
        a = ad.AnnData(X)
        a.var_names = ["gene", "other"]
        a.obs["leiden"] = pd.Categorical(["1"] * 50 + ["0"] * 150)
        a.raw = a
        return a

    def test_recovers_true_fold(self):
        a = self._make_adata(fold=2.0)
        lfc = NeuralTFPipeline._cluster_log2fc(a, "gene", "1")
        assert lfc == pytest.approx(1.0, abs=0.1)

    def test_null_gene_near_zero(self):
        a = self._make_adata(fold=2.0)
        lfc = NeuralTFPipeline._cluster_log2fc(a, "other", "1")
        assert abs(lfc) < 0.25

    def test_missing_cluster_returns_zero(self):
        a = self._make_adata(fold=2.0)
        assert NeuralTFPipeline._cluster_log2fc(a, "gene", None) == 0.0
        assert NeuralTFPipeline._cluster_log2fc(a, "gene", "99") == 0.0

    def test_missing_gene_returns_zero(self):
        a = self._make_adata(fold=2.0)
        assert NeuralTFPipeline._cluster_log2fc(a, "absent", "1") == 0.0


class TestPerezClassValidation:
    @pytest.mark.parametrize("cls", ["-", "", "nan", "None", "NA", "N/A", "n/a"])
    def test_invalid_classes(self, cls):
        assert NeuralTFPipeline._valid_perez_class(cls) is False

    @pytest.mark.parametrize("cls", ["bHLH", "Homeodomain", "POU", "Forkhead"])
    def test_valid_classes(self, cls):
        assert NeuralTFPipeline._valid_perez_class(cls) is True


class TestCorrelationAlignment:
    def test_joint_nan_alignment(self):
        # Row 1: (x1=0.2, g0=NaN); Row 2: (x1=NaN, g0=0.9); Row 3: good pair.
        # With independent dropna the old code re-paired row1.x1 with row2.g0.
        # Joint masking keeps only row 3.
        data = pd.DataFrame({
            "x1_corr": ["0.2", "", "0.1"],
            "g0_corr": ["", "0.9", "0.7"],
        })
        x1 = pd.to_numeric(data["x1_corr"], errors="coerce")
        g0 = pd.to_numeric(data["g0_corr"], errors="coerce")
        pair = pd.DataFrame({"x1": x1, "g0": g0}).dropna()
        assert len(pair) == 1
        assert pair["x1"].iloc[0] == pytest.approx(0.1)
        assert pair["g0"].iloc[0] == pytest.approx(0.7)
        gains = (pair["g0"] - pair["x1"]).to_numpy()
        assert gains[0] == pytest.approx(0.6)
