"""Tests for bioforge.plugins."""
from typing import Any

import pytest

from bioforge.core.exceptions import PluginError
from bioforge.plugins.base import PluginBase, PluginBaseImpl
from bioforge.plugins.manager import PluginManager


class _GoodPlugin(PluginBaseImpl):
    name = "good"
    version = "1.0.0"

    def initialize(self, config: Any) -> None:  # noqa: D401
        self.config = config


class _BadPlugin:
    """Missing 'name' and 'version' — fails PluginBase check."""
    def initialize(self, config: Any) -> None:
        pass


class _InitFailsPlugin(PluginBaseImpl):
    name = "initfails"
    version = "0.1.0"

    def initialize(self, config: Any) -> None:
        raise RuntimeError("init blew up")


class _ConstructorFailsPlugin(PluginBaseImpl):
    name = "ctorfails"
    version = "0.1.0"

    def __init__(self) -> None:
        raise RuntimeError("ctor failed")


class _DuplicateNamePlugin(PluginBaseImpl):
    name = "good"  # same name as _GoodPlugin
    version = "0.2.0"


def _fake_ep(value_cls, name: str):
    """Build a stub Entry-point-like object whose .load() returns a class."""
    class _Fake:
        def __init__(self, name: str, value: str) -> None:
            self.name = name
            self.value = value

        def load(self):
            return value_cls
    return _Fake(name, "x.y:Z")


def test_plugin_base_protocol_satisfied_by_impl() -> None:
    p = _GoodPlugin()
    assert isinstance(p, PluginBase)


def test_load_good_plugin() -> None:
    mgr = PluginManager()
    p = mgr.load(_fake_ep(_GoodPlugin, "good"), config={"x": 1})
    assert p.name == "good"
    assert p.version == "1.0.0"
    assert p.config == {"x": 1}
    assert mgr.names() == ["good"]
    assert len(mgr) == 1


def test_load_bad_plugin_missing_attrs() -> None:
    mgr = PluginManager()
    with pytest.raises(PluginError):
        mgr.load(_fake_ep(_BadPlugin, "bad"), config=None)


def test_init_failure_raises_plugin_error() -> None:
    mgr = PluginManager()
    with pytest.raises(PluginError):
        mgr.load(_fake_ep(_InitFailsPlugin, "initfails"), config=None)


def test_ctor_failure_raises_plugin_error() -> None:
    mgr = PluginManager()
    with pytest.raises(PluginError):
        mgr.load(_fake_ep(_ConstructorFailsPlugin, "ctorfails"), config=None)


def test_duplicate_name_raises() -> None:
    mgr = PluginManager()
    mgr.load(_fake_ep(_GoodPlugin, "good"), config=None)
    with pytest.raises(PluginError):
        mgr.load(_fake_ep(_DuplicateNamePlugin, "dup"), config=None)


def test_get_unknown_name_raises() -> None:
    mgr = PluginManager()
    with pytest.raises(PluginError):
        mgr.get("nope")


def test_iteration_yields_loaded_plugins() -> None:
    mgr = PluginManager()
    mgr.load(_fake_ep(_GoodPlugin, "good"), config=None)
    plugins = list(mgr)
    assert len(plugins) == 1
    assert plugins[0].name == "good"
