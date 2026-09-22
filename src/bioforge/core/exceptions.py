"""BioForge exception hierarchy.

All BioForge-specific errors inherit from :class:`BioForgeError`, allowing
callers to catch framework errors narrowly without masking unrelated bugs.
"""

from __future__ import annotations


class BioForgeError(Exception):
    """Base class for all BioForge-raised errors."""


# ----------------------------------------------------------------------------
# Configuration subsystem
# ----------------------------------------------------------------------------
class ConfigError(BioForgeError):
    """Raised when a configuration file cannot be loaded or is invalid."""


# ----------------------------------------------------------------------------
# Dataset subsystem
# ----------------------------------------------------------------------------
class DatasetError(BioForgeError):
    """Raised when a dataset cannot be located or is structurally invalid."""


# ----------------------------------------------------------------------------
# Plugin subsystem
# ----------------------------------------------------------------------------
class PluginError(BioForgeError):
    """Raised when a plugin fails to load, register, or validate."""


# ----------------------------------------------------------------------------
# Project subsystem
# ----------------------------------------------------------------------------
class ProjectError(BioForgeError):
    """Raised when a research project cannot be located or is misconfigured."""


__all__ = [
    "BioForgeError",
    "ConfigError",
    "DatasetError",
    "PluginError",
    "ProjectError",
]
