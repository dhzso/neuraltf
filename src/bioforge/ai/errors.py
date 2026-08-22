"""Exceptions raised by the AI Layer."""
from __future__ import annotations


class AIError(Exception):
    """Base class for AI-layer errors."""


class AIProviderNotConfiguredError(AIError):
    """No AI provider configured (no API key found)."""


class AIRequestError(AIError):
    """The underlying HTTP transport returned an error."""
    def __init__(self, message: str, status_code: int | None = None) -> None:
        super().__init__(message)
        self.status_code = status_code
