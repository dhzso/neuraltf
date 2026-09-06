"""BioForge research project manager.

Research projects live under ``<repo>/projects/``. Each project is a
subdirectory with a standardized layout (mirrors the NeuralTF scaffold
from the project charter)::

    projects/
        <project>/
            data/
            docs/
            figures/
            logs/
            notebooks/
            results/
            scripts/

Responsibilities:
- Validate that a project directory follows the standard layout.
- List available projects.
- Resolve a project by name and report its canonical paths.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bioforge.core.exceptions import ProjectError

# Canonical subdirectory structure every research project is expected to
# have. Used for both validation (do these directories exist?) and for
# creating new projects without sprawling ad-hoc structure.
PROJECT_SUBDIRS: tuple[str, ...] = (
    "data",
    "docs",
    "figures",
    "logs",
    "notebooks",
    "results",
    "scripts",
)


@dataclass(frozen=True)
class Project:
    """A resolved research project."""

    name: str
    root: Path

    @property
    def data(self) -> Path:
        return self.root / "data"

    @property
    def docs(self) -> Path:
        return self.root / "docs"

    @property
    def figures(self) -> Path:
        return self.root / "figures"

    @property
    def logs(self) -> Path:
        return self.root / "logs"

    @property
    def notebooks(self) -> Path:
        return self.root / "notebooks"

    @property
    def results(self) -> Path:
        return self.root / "results"

    @property
    def scripts(self) -> Path:
        return self.root / "scripts"


class ProjectManager:
    """Manage research projects under ``projects/``."""

    def __init__(self, root: str | Path):
        """Initialize the manager.

        Parameters
        ----------
        root
            Path to the BioForge repository root (the directory that
            contains ``projects/``).

        Raises
        ------
        ProjectError
            If ``<root>/projects`` does not exist.
        """
        self._root = Path(root).resolve()
        self._projects_dir = self._root / "projects"
        if not self._projects_dir.is_dir():
            raise ProjectError(
                f"Projects directory not found: {self._projects_dir}"
            )

    def list(self) -> list[str]:
        """Return the names of all immediate project subdirectories."""
        return sorted(
            p.name
            for p in self._projects_dir.iterdir()
            if p.is_dir() and not p.name.startswith(".")
        )

    def resolve(self, name: str) -> Project:
        """Resolve a project by name, validating its canonical layout.

        Raises :class:`ProjectError` if the project does not exist or is
        missing required subdirectories.
        """
        root = self._projects_dir / name
        if not root.is_dir():
            raise ProjectError(f"Project '{name}' not found at {root}.")
        missing = [s for s in PROJECT_SUBDIRS if not (root / s).is_dir()]
        if missing:
            raise ProjectError(
                f"Project '{name}' is missing standard subdirs: {missing}."
            )
        return Project(name=name, root=root)

    def create(self, name: str, exist_ok: bool = False) -> Project:
        """Create a new project scaffold under ``projects/``.

        Parameters
        ----------
        name
            Project directory name.
        exist_ok
            If ``True`` and the project already exists with the full
            standard layout, return it; otherwise raise.
        """
        root = self._projects_dir / name
        if root.exists():
            if not exist_ok:
                raise ProjectError(f"Project '{name}' already exists.")
            return self.resolve(name)
        for sub in PROJECT_SUBDIRS:
            (root / sub).mkdir(parents=True, exist_ok=True)
        return Project(name=name, root=root)

    def __repr__(self) -> str:  # pragma: no cover
        return f"ProjectManager(projects_dir={self._projects_dir!s})"


__all__ = ["PROJECT_SUBDIRS", "Project", "ProjectManager"]
