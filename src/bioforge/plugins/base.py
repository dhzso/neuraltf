"""BioForge plugin base protocol.

Defines the interface all BioForge plugins must implement. Plugins are
discovered via the ``bioforge.plugins`` entry-point group in
``pyproject.toml``::

    [project.entry-points."bioforge.plugins"]
    my_plugin = "my_package.my_module:MyPlugin"

A plugin instance is created once at registration time and receives no
constructor arguments — all configuration is read from the BioForge config
by the plugin itself.
"""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class PluginBase(Protocol):
    """Minimal plugin contract.

    Concrete plugins may subclass ``PluginBaseImpl`` (provided below) or
    simply implement these attributes directly. The runtime-checkable
    Protocol keeps the contract loose (duck-typed) while still allowing
    isinstance checks at manager time.

    Required attributes
    -------------------
    name : str
        Unique, short identifier for the plugin (used in logs and CLI).
    version : str
        Plugin version string (PEP 440 recommended).
    """

    name: str
    version: str

    def initialize(self, config: Any) -> None:
        """Called once after discovery, with the loaded BioForge config.

        Plugins may pull whatever they need from ``config`` and store it
        for later use. Raises should be propagated as
        :class:`bioforge.core.exceptions.PluginError`.
        """
        ...


class PluginBaseImpl:
    """Concrete base class plugins may extend for convenience.

    Provides default ``__init__``, name, and version attributes that
    subclasses override via class attributes. Use this when you want
    inheritance-based plugins rather than duck-typed Protocols.
    """

    name: str = "plugin"
    version: str = "0.0.0"

    def initialize(self, config: Any) -> None:  # noqa: D401
        """Default no-op initializer. Subclasses override."""
        return None


__all__ = ["PluginBase", "PluginBaseImpl"]
