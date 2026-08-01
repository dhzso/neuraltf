"""Smoke tests for bioforge.ui (Layer 10).

We test the pure-logic helpers under app.py directly so the tests don't
need a live Streamlit context.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from bioforge.ui.app import _pipeline_status, _runs_root, _repo_root


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
