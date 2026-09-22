"""BioForge logging setup.

Provides :func:`configure_logging`, which installs a single console handler
(and an optional file handler) on the root ``bioforge`` logger. The function
is idempotent: subsequent calls update the level/handlers without
duplicating handlers.

Typical usage::

    from bioforge.core.config import load_config
    from bioforge.core.logging import configure_logging

    cfg = load_config()
    configure_logging(cfg.logging)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from bioforge.core.config import LoggingConfig

# Single module-level logger instance for the package.
_LOGGER_NAME = "bioforge"


def _clear_handlers(logger: logging.Logger) -> None:
    """Remove all existing handlers from a logger."""
    for h in list(logger.handlers):
        logger.removeHandler(h)


def configure_logging(cfg: Optional[LoggingConfig] = None) -> logging.Logger:
    """Configure the ``bioforge`` logger.

    Parameters
    ----------
    cfg
        Logging configuration. If ``None``, defaults are used (INFO level,
        no file handler).

    Returns
    -------
    logging.Logger
        The configured ``bioforge`` logger.
    """
    if cfg is None:
        cfg = LoggingConfig()

    logger = logging.getLogger(_LOGGER_NAME)
    _clear_handlers(logger)
    logger.setLevel(cfg.level.upper())
    logger.propagate = False  # avoid duplicate output via root logger

    fmt = logging.Formatter(
        fmt="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler (always present)
    console = logging.StreamHandler()
    console.setFormatter(fmt)
    logger.addHandler(console)

    # Optional file handler
    if cfg.file:
        file_path = Path(cfg.file)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(file_path, encoding="utf-8")
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    return logger


def get_logger(name: str) -> logging.Logger:
    """Return a child logger of the ``bioforge`` package logger."""
    if not name:
        return logging.getLogger(_LOGGER_NAME)
    if name.startswith(_LOGGER_NAME + "."):
        return logging.getLogger(name)
    return logging.getLogger(f"{_LOGGER_NAME}.{name}")


__all__ = ["configure_logging", "get_logger"]
