"""Tests for bioforge.core.config."""
import os
from pathlib import Path

import pytest

from bioforge.core.config import (
    BioForgeConfig,
    DatasetsConfig,
    LoggingConfig,
    load_config,
)
from bioforge.core.exceptions import ConfigError


def _write_yaml(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def test_defaults_when_no_file(monkeypatch, tmp_path) -> None:
    monkeypatch.chdir(tmp_path)
    # No config file present, no env vars → defaults returned
    for key in list(os.environ):
        if key.startswith("BIOFORGE_"):
            monkeypatch.delenv(key, raising=False)
    cfg = load_config()
    assert isinstance(cfg, BioForgeConfig)
    assert cfg.project == "BioForge"
    assert isinstance(cfg.logging, LoggingConfig)
    assert cfg.logging.level == "INFO"
    assert isinstance(cfg.datasets, DatasetsConfig)
    assert cfg.datasets.root == "datasets"


def test_loads_yaml(monkeypatch, tmp_path) -> None:
    _write_yaml(
        tmp_path / "config" / "bioforge.yaml",
        "project: NeuralTF\n"
        "logging:\n  level: DEBUG\n  file: logs/bf.log\n"
        "datasets:\n  root: data\n",
    )
    monkeypatch.chdir(tmp_path)
    for key in list(os.environ):
        if key.startswith("BIOFORGE_"):
            monkeypatch.delenv(key, raising=False)
    cfg = load_config()
    assert cfg.project == "NeuralTF"
    assert cfg.logging.level == "DEBUG"
    assert cfg.logging.file == "logs/bf.log"
    assert cfg.datasets.root == "data"


def test_env_overrides_file(monkeypatch, tmp_path) -> None:
    _write_yaml(
        tmp_path / "config" / "bioforge.yaml",
        "project: NeuralTF\nlogging:\n  level: INFO\n",
    )
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("BIOFORGE_LOGGING_LEVEL", "WARNING")
    cfg = load_config()
    assert cfg.project == "NeuralTF"  # from file
    assert cfg.logging.level == "WARNING"  # overridden by env


def test_env_overrides_when_no_file(monkeypatch, tmp_path) -> None:
    monkeypatch.chdir(tmp_path)
    for key in list(os.environ):
        if key.startswith("BIOFORGE_"):
            monkeypatch.delenv(key, raising=False)
    monkeypatch.setenv("BIOFORGE_PROJECT", "EnvProject")
    cfg = load_config()
    assert cfg.project == "EnvProject"


def test_invalid_yaml_raises(monkeypatch, tmp_path) -> None:
    _write_yaml(
        tmp_path / "config" / "bioforge.yaml",
        "project: : bad: yaml: here\n",
    )
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ConfigError):
        load_config()


def test_non_mapping_root_raises(monkeypatch, tmp_path) -> None:
    _write_yaml(
        tmp_path / "config" / "bioforge.yaml",
        "- just\n- a\n- list\n",
    )
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ConfigError):
        load_config()


def test_config_is_immutable() -> None:
    cfg = BioForgeConfig()
    with pytest.raises(Exception):
        cfg.project = "Mutated"  # frozen dataclass
