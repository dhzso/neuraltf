"""BioForge configuration manager.

Loads configuration from a YAML file (default ``config/bioforge.yaml`` at
the repository root), then overlays environment variable overrides of the
form ``BIOFORGE_<SECTION>_<KEY>``. Configuration is returned as a frozen
dataclass tree for safe read-only access by callers.

Example YAML::

    project: NeuralTF
    datasets:
      root: datasets
    logging:
      level: INFO
      file: null
"""

from __future__ import annotations

import dataclasses
import os
from pathlib import Path
from typing import Any

import yaml

from bioforge.core.exceptions import ConfigError


# ----------------------------------------------------------------------------
# Configuration dataclasses
# ----------------------------------------------------------------------------
@dataclasses.dataclass(frozen=True)
class LoggingConfig:
    """Logging subsystem configuration."""

    level: str = "INFO"
    file: str | None = None


@dataclasses.dataclass(frozen=True)
class DatasetsConfig:
    """Dataset layout configuration (paths relative to repo root)."""

    root: str = "datasets"
    raw: str = "raw"
    processed: str = "processed"
    reference: str = "reference"
    cache: str = "cache"


@dataclasses.dataclass(frozen=True)
class BioForgeConfig:
    """Top level BioForge configuration."""

    project: str = "BioForge"
    logging: LoggingConfig = dataclasses.field(default_factory=LoggingConfig)
    datasets: DatasetsConfig = dataclasses.field(default_factory=DatasetsConfig)


# ----------------------------------------------------------------------------
# Loader
# ----------------------------------------------------------------------------
_DEFAULT_CONFIG = BioForgeConfig()


def _merge(target: dict[str, Any], overrides: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge overrides into target (mutates target)."""
    for k, v in overrides.items():
        if isinstance(v, dict) and isinstance(target.get(k), dict):
            _merge(target[k], v)
        else:
            target[k] = v
    return target


def _env_overrides() -> dict[str, Any]:
    """Build a nested dict from BIOFORGE_<SECTION>_<KEY> env vars.

    Top-level fields (for example ``BIOFORGE_PROJECT``) are placed directly
    in the returned dict by lower-casing the remainder after the
    ``BIOFORGE_`` prefix. Section-key overrides such as
    ``BIOFORGE_LOGGING_LEVEL`` are nested under their section.
    """
    out: dict[str, Any] = {}
    for name, value in os.environ.items():
        if not name.startswith("BIOFORGE_"):
            continue
        remainder = name[len("BIOFORGE_"):].lower()
        parts = remainder.split("_", 1)
        if len(parts) == 1:
            # Top-level field, e.g. BIOFORGE_PROJECT
            key = parts[0]
            # Only override if this is a known top-level field.
            if key in BioForgeConfig.__dataclass_fields__:
                out[key] = value
            continue
        if len(parts) == 2:
            section, key = parts
            # Only nest under known sections.
            if section in BioForgeConfig.__dataclass_fields__:
                out.setdefault(section, {})[key] = value
            continue
    return out


def _build_config(data: dict[str, Any]) -> BioForgeConfig:
    """Construct a BioForgeConfig from a flat dict."""
    raw = dict(data)
    # Nested sections
    logging_data = raw.pop("logging", {}) or {}
    datasets_data = raw.pop("datasets", {}) or {}
    plain = {
        k: v
        for k, v in raw.items()
        if k in BioForgeConfig.__dataclass_fields__
    }
    return BioForgeConfig(
        logging=LoggingConfig(**logging_data),
        datasets=DatasetsConfig(**datasets_data),
        **plain,
    )


def load_config(path: str | os.PathLike[str] | None = None) -> BioForgeConfig:
    """Load BioForge configuration from YAML, with environment overrides.

    Parameters
    ----------
    path
        Path to a YAML configuration file. If ``None``, the default
        ``config/bioforge.yaml`` (resolved relative to the current working
        directory) is used. If the file does not exist, defaults are returned
        with environment overrides applied on top.

    Returns
    -------
    BioForgeConfig
        Frozen, immutable configuration tree.

    Raises
    ------
    ConfigError
        If the YAML file exists but cannot be parsed or contains invalid
        structure.
    """
    data: dict[str, Any] = {}
    target = Path(path) if path else Path.cwd() / "config" / "bioforge.yaml"
    if target.is_file():
        try:
            with target.open("r", encoding="utf-8") as f:
                loaded = yaml.safe_load(f) or {}
        except yaml.YAMLError as exc:
            raise ConfigError(f"Failed to parse {target}: {exc}") from exc
        if not isinstance(loaded, dict):
            raise ConfigError(
                f"Configuration root must be a mapping; got {type(loaded).__name__}"
            )
        data = loaded
    # Apply environment overrides on top of file values
    _merge(data, _env_overrides())
    return _build_config(data)


__all__ = [
    "BioForgeConfig",
    "LoggingConfig",
    "DatasetsConfig",
    "load_config",
]
