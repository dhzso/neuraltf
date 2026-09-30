"""Tests for bioforge.projects.manager."""
from pathlib import Path

import pytest

from bioforge.core.exceptions import ProjectError
from bioforge.projects.manager import (
    PROJECT_SUBDIRS,
    Project,
    ProjectManager,
)


def _make_repo(root: Path) -> None:
    (root / "projects").mkdir(parents=True, exist_ok=True)


def test_subdirs_match_charter() -> None:
    assert PROJECT_SUBDIRS == (
        "data", "docs", "figures", "logs",
        "notebooks", "results", "scripts",
    )


def test_missing_projects_dir_raises(tmp_path) -> None:
    with pytest.raises(ProjectError):
        ProjectManager(tmp_path)


def test_list_empty(tmp_path) -> None:
    _make_repo(tmp_path)
    mgr = ProjectManager(tmp_path)
    assert mgr.list() == []


def test_list_returns_dirs_only(tmp_path) -> None:
    _make_repo(tmp_path)
    (tmp_path / "projects" / "Alpha").mkdir()
    (tmp_path / "projects" / "Beta").mkdir()
    (tmp_path / "projects" / ".hidden").mkdir()
    (tmp_path / "projects" / "file.txt").touch()
    mgr = ProjectManager(tmp_path)
    assert mgr.list() == ["Alpha", "Beta"]


def test_resolve_missing_project_raises(tmp_path) -> None:
    _make_repo(tmp_path)
    mgr = ProjectManager(tmp_path)
    with pytest.raises(ProjectError):
        mgr.resolve("nope")


def test_resolve_missing_subdir_raises(tmp_path) -> None:
    _make_repo(tmp_path)
    (tmp_path / "projects" / "Alpha").mkdir()
    mgr = ProjectManager(tmp_path)
    with pytest.raises(ProjectError):
        mgr.resolve("Alpha")  # missing all subdirs


def test_resolve_full_layout(tmp_path) -> None:
    _make_repo(tmp_path)
    (tmp_path / "projects" / "Alpha").mkdir()
    for sub in PROJECT_SUBDIRS:
        (tmp_path / "projects" / "Alpha" / sub).mkdir()
    mgr = ProjectManager(tmp_path)
    p = mgr.resolve("Alpha")
    assert isinstance(p, Project)
    assert p.name == "Alpha"
    assert p.data.is_dir()
    assert p.notebooks.is_dir()
    assert p.scripts.is_dir()


def test_create_scaffold(tmp_path) -> None:
    _make_repo(tmp_path)
    mgr = ProjectManager(tmp_path)
    p = mgr.create("Gamma")
    assert p.root.is_dir()
    assert p.data.is_dir()
    for sub in PROJECT_SUBDIRS:
        assert (p.root / sub).is_dir()


def test_create_existing_raises(tmp_path) -> None:
    _make_repo(tmp_path)
    mgr = ProjectManager(tmp_path)
    mgr.create("Alpha")
    with pytest.raises(ProjectError):
        mgr.create("Alpha")


def test_create_existing_with_exist_ok_returns_resolved(tmp_path) -> None:
    _make_repo(tmp_path)
    mgr = ProjectManager(tmp_path)
    mgr.create("Alpha")
    p = mgr.create("Alpha", exist_ok=True)
    assert p.name == "Alpha"
    assert p.data.is_dir()
