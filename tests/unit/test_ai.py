"""Unit tests for bioforge.ai (Layer 6)."""
from __future__ import annotations

import json

import pytest

from bioforge.ai import (
    AIAssistant,
    AIProviderNotConfiguredError,
    ChatMessage,
    ChatResponse,
    OpenAICompatClient,
    StubAssistant,
    build_assistant,
    list_tools,
    lookup_gene,
    register_tool,
    summarize_candidates,
)
from bioforge.ai.tools import get_tool


# ---------------------------------------------------------------------------
# build_assistant
# ---------------------------------------------------------------------------


def test_build_assistant_falls_back_to_stub_when_no_api_key(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("BIOFORGE_AI_API_KEY", raising=False)
    a = build_assistant()
    assert isinstance(a, StubAssistant)


def test_build_assistant_uses_openai_client_when_api_key_set(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("BIOFORGE_AI_API_KEY", "test-key")
    monkeypatch.setenv("BIOFORGE_AI_BASE_URL", "stub://provider")
    monkeypatch.setenv("BIOFORGE_AI_MODEL", "test-model")
    a = build_assistant()
    assert isinstance(a, OpenAICompatClient)
    assert a.api_key == "test-key"
    assert a.model == "test-model"


# ---------------------------------------------------------------------------
# StubAssistant
# ---------------------------------------------------------------------------


def test_stub_assistant_returns_canned_response() -> None:
    a = StubAssistant()
    r = a.complete([ChatMessage(role="user", content="hello world")])
    assert isinstance(r, ChatResponse)
    assert "bioforge-stub" in r.content
    assert "hello world" in r.content


def test_stub_assistant_require_live_raises() -> None:
    a = StubAssistant()
    with pytest.raises(AIProviderNotConfiguredError):
        a.complete([ChatMessage(role="user", content="x")], require_live=True)


# ---------------------------------------------------------------------------
# OpenAICompatClient — uses the `stub://` base scheme so no network
# ---------------------------------------------------------------------------


def test_openai_compat_stub_scheme_returns_canonical_response() -> None:
    a = OpenAICompatClient(base_url="stub://provider", api_key="k", model="m-1")
    r = a.complete([ChatMessage(role="user", content="hi")])
    assert r.model == "m-1"
    assert "stub-openai-compat" in r.content
    assert r.raw["x-roles"] == ["user"]


# ---------------------------------------------------------------------------
# AI tools registry
# ---------------------------------------------------------------------------


def test_tool_registry_round_trip() -> None:
    # Force-register a test-only tool.
    @register_tool("test_temp_tool")
    def fn():
        return "ok"
    assert "test_temp_tool" in list_tools()
    assert get_tool("test_temp_tool")() == "ok"


def test_lookup_gene_returns_no_match_when_bridge_missing(tmp_path) -> None:
    out = lookup_gene("soxB", bridge_path=str(tmp_path / "missing.csv"))
    payload = json.loads(out)
    assert payload["ok"] is False


def test_lookup_gene_matches_by_name(tmp_path) -> None:
    import pandas as pd
    df = pd.DataFrame({
        "gene_name": ["soxB", "myoD"],
        "v6_id": ["v6_1", "v6_2"],
        "v4_id": ["v4_1", "v4_2"],
    })
    p = tmp_path / "bridge.csv"
    df.to_csv(p, index=False)
    out = lookup_gene("soxB", bridge_path=str(p))
    payload = json.loads(out)
    assert payload["ok"] is True
    assert payload["gene"]["v6_id"] == "v6_1"


def test_lookup_gene_matches_by_v6_id(tmp_path) -> None:
    import pandas as pd
    df = pd.DataFrame({
        "gene_name": ["soxB"],
        "v6_id": ["v6_1"],
        "v4_id": ["v4_1"],
    })
    p = tmp_path / "bridge.csv"
    df.to_csv(p, index=False)
    # v6_id lookup uses column-name equality; both column header and row match the query.
    out = lookup_gene("v6_1", bridge_path=str(p))
    payload = json.loads(out)
    assert payload["ok"] is True


def test_summarize_candidates_reads_existing_csv(tmp_path) -> None:
    p = tmp_path / "rank.csv"
    p.write_text(
        "gene_id,gene_name,tier,integrated_score\n"
        "v6_1,soxB,high,0.9\n"
        "v6_2,myoD,medium,0.4\n"
        "v6_3,unknown,low,0.1\n"
    )
    out = summarize_candidates(str(p), top_n=2)
    payload = json.loads(out)
    assert payload["ok"] is True
    assert payload["n_total"] == 3
    assert len(payload["top"]) == 2
    assert payload["top"][0]["gene_id"] == "v6_1"


def test_summarize_candidates_file_missing(tmp_path) -> None:
    out = summarize_candidates(str(tmp_path / "nope.csv"))
    payload = json.loads(out)
    assert payload["ok"] is False


# ---------------------------------------------------------------------------
# Subclassing AIAssistant should be straightforward (smoke test)
# ---------------------------------------------------------------------------


def test_assistant_subclass_complete_contract() -> None:
    class Dummy(AIAssistant):
        name = "dummy"
        def complete(self, messages):
            return ChatResponse(content="ok", model="dummy", elapsed_s=0.0)

    r = Dummy().complete([ChatMessage(role="user", content="x")])
    assert r.content == "ok"
