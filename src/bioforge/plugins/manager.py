"""BioForge plugin manager.

Discovers, loads, and validates plugins registered under the
``bioforge.plugins`` entry-point group. Plugin instances are created once
per manager and indexed by ``name``.

Usage::

    from bioforge.plugins import PluginManager

    mgr = PluginManager()
    mgr.load_all(config)
    for plugin in mgr:
        print(plugin.name, plugin.version)
"""

from __future__ import annotations

from importlib import metadata
from typing import Any, Iterator

from bioforge.core.exceptions import PluginError
from bioforge.plugins.base import PluginBase


_ENTRY_POINT_GROUP = "bioforge.plugins"


class PluginManager:
    """Registry of discovered and initialized BioForge plugins."""

    def __init__(self) -> None:
        self._plugins: dict[str, PluginBase] = {}

    # -- Discovery ---------------------------------------------------------
    def discover(self) -> list[metadata.EntryPoint]:
        """Return the list of entry points in the ``bioforge.plugins`` group.

        This only enumerates entry points; it does not load the plugin code.
        """
        return list(
            metadata.entry_points(group=_ENTRY_POINT_GROUP)
        )

    # -- Loading -----------------------------------------------------------
    def load(self, entry_point: metadata.EntryPoint, config: Any) -> PluginBase:
        """Load a single plugin from an entry point and initialize it.

        The plugin class is instantiated with no arguments, then ``initialize``
        is called with the supplied config. Failures raise
        :class:`PluginError`.
        """
        try:
            plugin_cls = entry_point.load()
        except Exception as exc:  # noqa: BLE001 - surface as PluginError
            raise PluginError(
                f"Failed to load plugin '{entry_point.name}': {exc}"
            ) from exc
        try:
            plugin = plugin_cls()
        except Exception as exc:  # noqa: BLE001
            raise PluginError(
                f"Plugin '{entry_point.name}' constructor failed: {exc}"
            ) from exc
        if not hasattr(plugin, "name") or not hasattr(plugin, "version"):
            raise PluginError(
                f"Plugin '{entry_point.name}' does not satisfy PluginBase "
                f"(missing 'name' or 'version' attribute)."
            )
        try:
            plugin.initialize(config)
        except PluginError:
            raise
        except Exception as exc:  # noqa: BLE001
            raise PluginError(
                f"Plugin '{entry_point.name}' initialize() failed: {exc}"
            ) from exc
        if plugin.name in self._plugins:
            raise PluginError(
                f"Duplicate plugin name '{plugin.name}' "
                f"(from {entry_point.name})"
            )
        self._plugins[plugin.name] = plugin
        return plugin

    def load_all(self, config: Any) -> list[PluginBase]:
        """Discover and load every registered plugin.

        Returns the list of successfully loaded plugins. Plugins that fail
        to load are skipped (their error is recorded as a warning); a future
        hard-fail mode may be added. This soft-fail behavior keeps one bad
        plugin from breaking BioForge's startup.
        """
        loaded: list[PluginBase] = []
        for ep in self.discover():
            try:
                loaded.append(self.load(ep, config))
            except PluginError as exc:
                # Soft-fail: log via stderr-only logger so we never crash the host.
                import sys
                print(
                    f"[bioforge.plugins] WARNING: skipping plugin "
                    f"'{ep.name}': {exc}",
                    file=sys.stderr,
                )
        return loaded

    # -- Querying ----------------------------------------------------------
    def get(self, name: str) -> PluginBase:
        """Return the plugin with the given unique name."""
        try:
            return self._plugins[name]
        except KeyError as exc:
            raise PluginError(f"No plugin named '{name}'.") from exc

    def names(self) -> list[str]:
        return sorted(self._plugins.keys())

    def __iter__(self) -> Iterator[PluginBase]:
        return iter(self._plugins.values())

    def __len__(self) -> int:
        return len(self._plugins)


__all__ = ["PluginManager"]
