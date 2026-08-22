"""Tests for bioforge.core.datasets."""
from pathlib import Path

import pytest

from bioforge.core.config import DatasetsConfig
from bioforge.core.datasets import (
    CATEGORIES,
    DatasetManager,
    DatasetPath,
)
from bioforge.core.exceptions import DatasetError


def _make_tree(root: Path) -> None:
    base = root / "datasets"
    base.mkdir(parents=True, exist_ok=True)
    for cat in CATEGORIES:
        (base / cat).mkdir(exist_ok=True)
    (base / "raw" / "AtlasA").mkdir(exist_ok=True)
    (base / "raw" / "AtlasA" / "data.h5ad").touch()


def test_categories_is_canonical_set() -> None:
    assert CATEGORIES == ("raw", "processed", "reference", "cache")


def test_manager_initialization(tmp_path) -> None:
    _make_tree(tmp_path)
    mgr = DatasetManager(tmp_path)
    assert mgr.base == (tmp_path / "datasets").resolve()


def test_missing_root_raises(tmp_path) -> None:
    with pytest.raises(DatasetError):
        DatasetManager(tmp_path / "nonexistent")


def test_unknown_category_raises(tmp_path) -> None:
    _make_tree(tmp_path)
    mgr = DatasetManager(tmp_path)
    with pytest.raises(DatasetError):
        mgr.category_path("bad")


def test_resolve_directory(tmp_path) -> None:
    _make_tree(tmp_path)
    mgr = DatasetManager(tmp_path)
    dp = mgr.resolve("raw", "AtlasA")
    assert isinstance(dp, DatasetPath)
    assert dp.category == "raw"
    assert dp.path.exists()


def test_resolve_nested_file(tmp_path) -> None:
    _make_tree(tmp_path)
    mgr = DatasetManager(tmp_path)
    dp = mgr.resolve("raw", "AtlasA/data.h5ad")
    assert dp.path.exists()
    assert dp.path.is_file()


def test_resolve_blocks_path_traversal(tmp_path) -> None:
    _make_tree(tmp_path)
    mgr = DatasetManager(tmp_path)
    with pytest.raises(DatasetError):
        mgr.resolve("raw", "../../etc/passwd")


def test_list_returns_sorted_children(tmp_path) -> None:
    _make_tree(tmp_path)
    (tmp_path / "datasets" / "raw" / "AtlasB").mkdir()
    (tmp_path / "datasets" / "raw" / ".hidden").mkdir()
    mgr = DatasetManager(tmp_path)
    items = mgr.list("raw")
    assert items == ["AtlasA", "AtlasB"]  # sorted, hidden excluded


def test_list_empty_category(tmp_path) -> None:
    _make_tree(tmp_path)
    mgr = DatasetManager(tmp_path)
    assert mgr.list("processed") == []


def test_config_relpaths(tmp_path) -> None:
    # Custom layout where 'raw' lives elsewhere (structurally)
    (tmp_path / "tsdata").mkdir()
    (tmp_path / "tsdata" / "rawdata").mkdir()
    cfg = DatasetsConfig(root="tsdata", raw="rawdata")
    mgr = DatasetManager(tmp_path, cfg)
    assert mgr.base.name == "tsdata"
    assert mgr.category_path("raw").name == "rawdata"
