"""Integration tests against the real mounted BioForge repository.

These tests run against the bind-mounted ``/workspace`` (the actual repo)
so they validate integration with the real directory layout rather than
synthetic temp fixtures. They require the container's working directory
to be ``/workspace``.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from bioforge.core.config import load_config
from bioforge.core.datasets import DatasetManager
from bioforge.projects.manager import ProjectManager


WORKSPACE = Path(os.environ.get("BIOFORGE_WORKSPACE", "/workspace"))


def _in_repo() -> bool:
    return (WORKSPACE / "pyproject.toml").is_file() and (WORKSPACE / "projects").is_dir()


pytestmark = pytest.mark.skipif(
    not _in_repo(),
    reason="Tests must run inside the BioForge container with /workspace mounted",
)


def test_repo_root_is_workspace() -> None:
    assert WORKSPACE.is_dir()
    assert (WORKSPACE / "pyproject.toml").is_file()


def test_real_datasets_tree_exists() -> None:
    base = WORKSPACE / "datasets"
    assert base.is_dir()
    for cat in ("raw", "processed", "reference", "cache"):
        assert (base / cat).is_dir()


def test_real_datasets_list_raw() -> None:
    mgr = DatasetManager(WORKSPACE)
    raw = mgr.list("raw")
    # Real datasets we expect from the supplied data
    assert any("GSE103633" in n for n in raw)
    assert any("GSE111764" in n for n in raw)


def test_real_datasets_reference_has_papers() -> None:
    mgr = DatasetManager(WORKSPACE)
    reference_list = mgr.list("reference")
    assert any("Papers" in n or "paper" in n.lower() for n in reference_list)


def test_real_neuraltf_project_present() -> None:
    mgr = ProjectManager(WORKSPACE)
    names = mgr.list()
    assert "NeuralTF" in names
    p = mgr.resolve("NeuralTF")
    # NeuralTF was created by the charter scaffold; all standard subdirs should be present
    assert p.data.is_dir()
    assert p.notebooks.is_dir()
    assert p.scripts.is_dir()
    assert p.results.is_dir()


def test_load_config_inside_repo_default_path() -> None:
    # The repo currently has no config/bioforge.yaml, so defaults are returned
    cfg = load_config()
    assert cfg.project == "BioForge"
    assert cfg.datasets.root == "datasets"


def test_cli_info_invokes_against_real_repo(monkeypatch) -> None:
    from click.testing import CliRunner

    from bioforge.cli.main import cli

    monkeypatch.chdir(WORKSPACE)
    runner = CliRunner()
    result = runner.invoke(cli, ["info"])
    assert result.exit_code == 0
    assert "bioforge" in result.output


def test_cli_datasets_list_against_real_repo(monkeypatch) -> None:
    from click.testing import CliRunner

    from bioforge.cli.main import cli

    monkeypatch.chdir(WORKSPACE)
    runner = CliRunner()
    result = runner.invoke(cli, ["datasets", "list", "--category", "raw"])
    assert result.exit_code == 0
    assert "GSE103633" in result.output or "GSE111764" in result.output


def test_cli_projects_list_against_real_repo(monkeypatch) -> None:
    from click.testing import CliRunner

    from bioforge.cli.main import cli

    monkeypatch.chdir(WORKSPACE)
    runner = CliRunner()
    result = runner.invoke(cli, ["projects", "list"])
    assert result.exit_code == 0
    assert "NeuralTF" in result.output
