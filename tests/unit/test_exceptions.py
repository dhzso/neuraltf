"""Tests for bioforge.core.exceptions."""
from bioforge.core.exceptions import (
    BioForgeError,
    ConfigError,
    DatasetError,
    PluginError,
    ProjectError,
)


def test_base_class_is_exception() -> None:
    assert issubclass(BioForgeError, Exception)


def test_subclasses_inherit_from_base() -> None:
    for exc in (ConfigError, DatasetError, PluginError, ProjectError):
        assert issubclass(exc, BioForgeError)


def test_subclasses_are_distinct() -> None:
    assert ConfigError is not DatasetError
    assert PluginError is not ProjectError


def test_can_raise_and_catch() -> None:
    try:
        raise ConfigError("boom")
    except BioForgeError as exc:
        assert str(exc) == "boom"
