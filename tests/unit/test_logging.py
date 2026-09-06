"""Tests for bioforge.core.logging."""
import logging
from pathlib import Path

from bioforge.core.config import LoggingConfig
from bioforge.core.logging import configure_logging, get_logger


def test_configure_logging_returns_logger() -> None:
    logger = configure_logging(LoggingConfig(level="INFO"))
    assert logger.name == "bioforge"
    assert logger.level == logging.INFO


def test_reconfigure_does_not_duplicate_handlers() -> None:
    configure_logging(LoggingConfig(level="INFO"))
    configure_logging(LoggingConfig(level="DEBUG"))
    logger = logging.getLogger("bioforge")
    assert len(logger.handlers) == 1


def test_file_handler(tmp_path) -> None:
    log_file = tmp_path / "subdir" / "bf.log"
    cfg = LoggingConfig(level="INFO", file=str(log_file))
    logger = configure_logging(cfg)
    test_logger = get_logger("test")
    test_logger.info("hello-xyz")
    # Flush handlers
    for h in logger.handlers:
        h.flush()
    assert log_file.is_file()
    content = log_file.read_text(encoding="utf-8")
    assert "hello-xyz" in content


def test_get_logger_child() -> None:
    name = "subsystem.foo"
    lg = get_logger(name)
    assert lg.name == f"bioforge.{name}"
    # Already prefixed is not re-prefixed
    lg2 = get_logger("bioforge.bar")
    assert lg2.name == "bioforge.bar"
