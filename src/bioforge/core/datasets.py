"""BioForge dataset manager.

Resolves dataset paths within the BioForge on-disk layout::

    <repo_root>/
        datasets/
        ├── raw/       — immutable original datasets
        ├── processed/ — derived artifacts (AnnData, etc.)
        ├── reference/ — reference genomes, annotations, papers
        └── cache/     — transient cache (gitignored)

The :class:`DatasetManager`` knows about this layout, validates that
directories exist, and offers :meth:`resolve` for accessing a named
dataset within any of the four categories.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bioforge.core.config import DatasetsConfig
from bioforge.core.exceptions import DatasetError

# Canonical category names. Locking these to a constant set keeps logs and
# downstream code consistent — no string-typed freelancing.
CATEGORIES: tuple[str, ...] = ("raw", "processed", "reference", "cache")


@dataclass(frozen=True)
class DatasetPath:
    """A resolved dataset location."""

    category: str
    name: str
    path: Path

    def exists(self) -> bool:
        return self.path.exists()


class DatasetManager:
    """Resolve and validate dataset paths under the BioForge layout."""

    def __init__(self, root: str | Path, config: DatasetsConfig | None = None):
        """Initialize the manager.

        Parameters
        ----------
        root
            Path to the repository root (the directory containing the
            ``datasets/`` tree).
        config
            Dataset layout configuration. If ``None``, defaults are used.

        Raises
        ------
        DatasetError
            If the ``<root>/<config.root>`` directory does not exist.
        """
        self._root = Path(root).resolve()
        self._config = config or DatasetsConfig()
        self._base = self._root / self._config.root
        if not self._base.is_dir():
            raise DatasetError(
                f"Dataset root not found: {self._base}. "
                f"Did you mean to create it?"
            )

    # -- Path accessors ----------------------------------------------------
    @property
    def base(self) -> Path:
        """Root of the datasets tree (e.g. ``<repo>/datasets``)."""
        return self._base

    def category_path(self, category: str) -> Path:
        """Return the path for one of the four canonical categories."""
        if category not in CATEGORIES:
            raise DatasetError(
                f"Unknown dataset category '{category}'. "
                f"Expected one of {CATEGORIES}."
            )
        # The .root attribute holds the top-level dir; sub-dirs use private
        # attribute names so we map them defensively.
        sub = getattr(self._config, category)
        return self._base / sub

    def resolve(self, category: str, name: str) -> DatasetPath:
        """Resolve a named dataset inside a category.

        ``name`` may be either a directory name (e.g. ``"GSE103633_GEO_Plass_atlas"``)
        or a relative path (e.g. ``"GSE103633_GEO_Plass_atlas/series_matrix.txt.gz"``).
        The returned :class:`DatasetPath` is guaranteed to point inside the
        requested category but is NOT required to exist on disk (the caller
        decides whether existence is mandatory).
        """
        cat_root = self.category_path(category)
        candidate = (cat_root / name).resolve()
        # Guard against path traversal: candidate must start with cat_root.
        try:
            candidate.relative_to(cat_root)
        except ValueError as exc:
            raise DatasetError(
                f"Resolved path '{candidate}' escapes category '{category}' "
                f"root '{cat_root}'."
            ) from exc
        return DatasetPath(category=category, name=name, path=candidate)

    # -- Discovery ---------------------------------------------------------
    def list(self, category: str) -> list[str]:
        """List immediate children of a category directory."""
        cat_root = self.category_path(category)
        if not cat_root.is_dir():
            return []
        return sorted(p.name for p in cat_root.iterdir() if not p.name.startswith("."))

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"DatasetManager(root={self._root!s}, base={self._base!s})"


__all__ = ["CATEGORIES", "DatasetPath", "DatasetManager"]
