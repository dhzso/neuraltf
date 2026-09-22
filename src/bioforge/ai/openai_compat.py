"""OpenAI-compatible HTTP client.

A single small client posts to ``${base_url}/chat/completions`` with
``Authorization: Bearer ${api_key}``. Works with OpenAI, OpenRouter,
Together, Groq, Anyscale, vLLM, Ollama's OpenAI shim, LM Studio, and
any other conforming provider.

We use only Python stdlib (``urllib.request`` + ``json``) so the AI Layer
adds zero runtime dependencies.
"""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from dataclasses import dataclass

from bioforge.ai.assistant import AIAssistant, ChatMessage, ChatResponse
from bioforge.ai.errors import AIRequestError
from bioforge.core.logging import get_logger

logger = get_logger("ai.openai_compat")


@dataclass
class OpenAICompatClient(AIAssistant):
    base_url: str
    api_key: str
    model: str
    timeout_s: float = 60.0
    name: str = "openai_compat"

    def complete(self, messages: list[ChatMessage]) -> ChatResponse:
        # Treat stub "OpenAI-compat" base URLs that begin with "stub://" as
        # test fixtures: instead of hitting the network, echo back a stable
        # canned string so unit tests stay hermetic.
        if str(self.base_url).startswith("stub://"):
            return ChatResponse(
                content=f"[stub-openai-compat] model={self.model}",
                model=self.model,
                elapsed_s=0.0,
                raw={"x-roles": [m.role for m in messages]},
            )

        url = f"{self.base_url.rstrip('/')}/chat/completions"
        payload = {
            "model": self.model,
            "messages": [
                {"role": m.role, "content": m.content} for m in messages
            ],
            "temperature": 0.2,
        }
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
            "Accept": "application/json",
        }
        body = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(url, data=body, headers=headers, method="POST")
        t0 = time.perf_counter()
        try:
            with urllib.request.urlopen(req, timeout=self.timeout_s) as resp:
                raw = json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            message = f"HTTP {exc.code} from {url}: {exc.read().decode('utf-8', 'replace')}"
            logger.error(message)
            raise AIRequestError(message, status_code=exc.code) from exc
        except urllib.error.URLError as exc:
            raise AIRequestError(f"network error: {exc.reason}") from exc
        elapsed = time.perf_counter() - t0
        try:
            content = raw["choices"][0]["message"]["content"]
        except (KeyError, IndexError, TypeError) as exc:
            raise AIRequestError(
                f"malformed response from {url}: {raw!r}"
            ) from exc
        logger.info("openai-compat completed in %.4fs (model=%s)", elapsed, self.model)
        return ChatResponse(content=content, model=raw.get("model", self.model),
                            elapsed_s=elapsed, raw=raw)
