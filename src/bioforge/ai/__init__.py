"""BioForge AI Layer (Layer 6).

A thin, provider-agnostic AI assistant interface. Concrete implementation
is OpenAI-compatible HTTP (works with OpenAI, OpenRouter, Together, Groq,
vLLM, Ollama's OpenAI shim, LM Studio, etc.). A built-in stub fallback lets
the rest of BioForge run with zero AI configuration, which matters for the
thesis reproducibility guarantee and for offline dev.

Consumers only ever call :meth:`AIAssistant.complete` and receive a
:class:`ChatResponse`. Configuring which assistant to use is done by
:func:`build_assistant` reading env vars and optional config files.
"""
from bioforge.ai.assistant import (
    AIAssistant,
    ChatMessage,
    ChatResponse,
    StubAssistant,
)
from bioforge.ai.errors import AIProviderNotConfiguredError
from bioforge.ai.openai_compat import OpenAICompatClient
from bioforge.ai.tools import list_tools, lookup_gene, register_tool, summarize_candidates

__all__ = [
    "AIAssistant",
    "ChatMessage",
    "ChatResponse",
    "StubAssistant",
    "OpenAICompatClient",
    "AIProviderNotConfiguredError",
    "build_assistant",
    "lookup_gene",
    "register_tool",
    "summarize_candidates",
    "list_tools",
]


def build_assistant(config: "dict | None" = None) -> AIAssistant:
    """Return the AI assistant the framework should use by default.

    Reads env vars (``BIOFORGE_AI_API_KEY`` et al.) first, then the
    ``ai`` section of the optional dict passed in. Falls back to
    :class:`StubAssistant` when no key is found.
    """
    import os

    cfg = {
        "base_url": os.environ.get(
            "BIOFORGE_AI_BASE_URL", "https://api.openai.com/v1"
        ),
        "api_key": os.environ.get("BIOFORGE_AI_API_KEY"),
        "model": os.environ.get("BIOFORGE_AI_MODEL", "gpt-4o-mini"),
    }
    if config:
        for k, v in config.get("ai", {}).items():
            cfg[k] = v
    if not cfg["api_key"]:
        return StubAssistant()
    return OpenAICompatClient(
        base_url=cfg["base_url"],
        api_key=cfg["api_key"],
        model=cfg["model"],
    )
