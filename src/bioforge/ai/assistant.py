"""Assistant interface and the Stub fallback.

The abstract :class:`AIAssistant` is intentionally minimal so any provider
that wants to be BioForge-compatible only needs to implement ``complete``.
The :class:`StubAssistant` returns deterministic canned responses, which
keeps unit tests hermetic and lets BioForge work offline without any AI
configuration while everything still dereferences.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional

from bioforge.ai.errors import AIProviderNotConfiguredError
from bioforge.core.logging import get_logger

logger = get_logger("ai.assistant")


@dataclass
class ChatMessage:
    role: str  # "system" | "user" | "assistant" | "tool"
    content: str
    name: Optional[str] = None


@dataclass
class ChatResponse:
    content: str
    model: str
    elapsed_s: float
    raw: dict = field(default_factory=dict)


class AIAssistant:
    """Abstract AI assistant interface."""

    name: str = "abstract"

    def complete(self, messages: list[ChatMessage]) -> ChatResponse:  # pragma: no cover
        raise NotImplementedError


class StubAssistant(AIAssistant):
    """Deterministic stub used when no AI provider is configured.

    Never raises (unless ``require_live=True`` is passed to ``complete``);
    returns a deterministic 'stubbed' reply so workflows degrade gracefully
    when AI is unavailable.
    """

    name: str = "stub"

    def complete(
        self,
        messages: list[ChatMessage],
        *,
        require_live: bool = False,
    ) -> ChatResponse:
        if require_live:
            raise AIProviderNotConfiguredError(
                "AI provider not configured; BIOFORGE_AI_API_KEY is unset."
            )
        t0 = time.perf_counter()
        # Echo a deterministic, intentional signal so tests can detect the stub.
        user_msg = next(
            (m for m in messages if m.role == "user"), None
        )
        body = "[bioforge-stub] AI provider not configured. "
        if user_msg is not None:
            snippet = user_msg.content[:140].replace("\n", " ")
            body += f"Would have answered: {snippet!r}"
        elapsed = time.perf_counter() - t0
        logger.debug("stub assistant replied (%.4fs)", elapsed)
        return ChatResponse(content=body, model="stub", elapsed_s=elapsed, raw={})
