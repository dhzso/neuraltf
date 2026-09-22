"""Smoke tests for bioforge.ui (Layer 10).

We test the pure-logic helpers under app.py directly so the tests don't
need a live Streamlit context.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from bioforge.ui.app import (
    _data_rebuild_readiness,
    _downstream_readiness,
    _pipeline_status,
    _post_pipeline_steps,
    _rebuild_data,
    _runs_root,
    _repo_root,
)


def test_pipeline_status_returns_dict(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # Build a fake repo root with the files needed
    (tmp_path / "datasets" / "processed").mkdir(parents=True)
    (tmp_path / "datasets" / "raw").mkdir(parents=True)
    (tmp_path / "projects" / "NeuralTF" / "data").mkdir(parents=True)
    (tmp_path / "datasets" / "processed" / "fincher_subsample.h5ad").write_bytes(b"")
    (tmp_path / "datasets" / "processed" / "plass_v6.h5ad").write_bytes(b"")

    monkeypatch.setenv("BIOFORGE_REPO_ROOT", str(tmp_path))
    status = _pipeline_status()
    assert isinstance(status, dict)
    assert status["fincher_h5ad"] is True
    assert status["plass_h5ad"] is True


def test_app_main_module_imports_without_streamlit_run() -> None:
    # Plain import sanity check; the deferred `_st()` indirection keeps the
    # module importable even if streamlit isn't installed.
    import bioforge.ui.app as app
    assert callable(app.main)
    assert callable(app.render_run_page)
    assert callable(app.render_results_page)
    assert callable(app.render_assistant_page)


def test_repo_root_respects_env_var(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("BIOFORGE_REPO_ROOT", str(tmp_path))
    assert _repo_root() == tmp_path


def test_runs_root_under_repo_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("BIOFORGE_REPO_ROOT", str(tmp_path))
    assert _runs_root() == tmp_path / "projects" / "NeuralTF" / "runs"


def test_downstream_readiness_gates_on_inputs(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("BIOFORGE_REPO_ROOT", str(tmp_path))
    ds = _downstream_readiness(tmp_path)
    assert ds["figures"] is True
    assert ds["prioritization"] is False
    assert ds["go_figures"] is False

    king = tmp_path / "datasets" / "raw" / "Supplementary_Data_ King_2024"
    king.mkdir(parents=True)
    (king / "1-s2.0-S2211124724001712-mmc4.xlsx").write_bytes(b"")
    (king / "1-s2.0-S2211124724001712-mmc5.xlsx").write_bytes(b"")
    (tmp_path / "datasets" / "processed").mkdir(parents=True)
    (tmp_path / "datasets" / "processed" / "planmine_annotations.parquet").write_bytes(b"")
    (tmp_path / "datasets" / "raw" / "go.obo").write_bytes(b"")

    ds = _downstream_readiness(tmp_path)
    assert ds["prioritization"] is True
    assert ds["go_figures"] is True


def test_post_pipeline_steps_gates_on_missing_inputs(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("BIOFORGE_REPO_ROOT", str(tmp_path))
    run_dir = tmp_path / "projects" / "NeuralTF" / "runs" / "x"
    run_dir.mkdir(parents=True)

    class _St:
        def __init__(self) -> None:
            self.captions: list[str] = []
            self.warnings: list[str] = []

        def caption(self, s: str) -> None:
            self.captions.append(s)

        def warning(self, s: str) -> None:
            self.warnings.append(s)

        def markdown(self, s: str) -> None:
            pass

    stub = _St()
    monkeypatch.setattr("bioforge.ui.app._st", lambda: stub)
    # No subprocess may run: every step must be gated off here.
    def _boom(*args, **kwargs):  # pragma: no cover - must not be called
        raise AssertionError("subprocess must not run when inputs are missing")
    monkeypatch.setattr("subprocess.run", _boom)

    box = type("Box", (), {"markdown": lambda self, s: None})()
    _post_pipeline_steps(tmp_path, run_dir, [], box)
    assert len(stub.captions) == 3  # figures, prioritization, GO figures all skipped


def test_data_rebuild_readiness_gates_on_raw_inputs(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("BIOFORGE_REPO_ROOT", str(tmp_path))
    assert _data_rebuild_readiness(tmp_path) == {
        "fincher": False, "plass": False, "bridge": False,
        "king_atlas": False, "planmine": False,
    }

    raw = tmp_path / "datasets" / "raw"
    king = raw / "Supplementary_Data_ King_2024"
    (raw / "GSE111764_GEO_Fincher_atlas").mkdir(parents=True)
    (raw / "GSE111764_GEO_Fincher_atlas"
     / "GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz").write_bytes(b"")
    (raw / "GSE103633_GEO_Plass_atlas").mkdir(parents=True)
    (raw / "GSE103633_GEO_Plass_atlas" / "GSE103633_RAW.tar").write_bytes(b"")
    (raw / "smed_20140614.mapping.rosettastone.2020.txt").write_bytes(b"")
    king.mkdir(parents=True)
    (king / "1-s2.0-S2211124724001712-mmc4.xlsx").write_bytes(b"")
    (king / "1-s2.0-S2211124724001712-mmc7.xlsx").write_bytes(b"")
    (tmp_path / "projects" / "NeuralTF" / "runs" / "pipeline_run").mkdir(parents=True)
    (tmp_path / "projects" / "NeuralTF" / "runs" / "pipeline_run"
     / "rank_neural.csv").write_text("gene_id\ng\na\n")

    assert _data_rebuild_readiness(tmp_path) == {
        "fincher": True, "plass": True, "bridge": True,
        "king_atlas": True, "planmine": True,
    }


def test_rebuild_data_gates_on_missing_inputs(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("BIOFORGE_REPO_ROOT", str(tmp_path))

    class _St:
        def __init__(self) -> None:
            self.captions: list[str] = []

        def caption(self, s: str) -> None:
            self.captions.append(s)

        def warning(self, s: str) -> None:
            pass

        def markdown(self, s: str) -> None:
            pass

    stub = _St()
    monkeypatch.setattr("bioforge.ui.app._st", lambda: stub)

    def _boom(*args, **kwargs):  # pragma: no cover - must not be called
        raise AssertionError("subprocess must not run when raw inputs are missing")
    monkeypatch.setattr("subprocess.run", _boom)

    box = type("Box", (), {"markdown": lambda self, s: None})()
    _rebuild_data(tmp_path, [], box)
    # all five rebuild steps skipped (fincher, plass, bridge, king atlas, planmine)
    assert len(stub.captions) == 5
