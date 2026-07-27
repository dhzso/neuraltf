"""Deterministic tools the AI assistant can call.

These are simple Python callables registered with the AI module. They give
an LLM-driven assistant agency to look up real BioForge data (a gene in a
bridge table, a ranked candidates CSV) without the LLM having direct I/O
access. A real tool-calling framework (function calling in OpenAI-compatible
providers) wraps these callables and returns their JSON output to the model.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Callable, Optional

from bioforge.core.logging import get_logger

logger = get_logger("ai.tools")


_TOOL_REGISTRY: dict[str, Callable[..., Any]] = {}


def register_tool(name: str) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Decorator: register a Python callable as a BioForge AI tool."""

    def deco(fn: Callable[..., Any]) -> Callable[..., Any]:
        if name in _TOOL_REGISTRY:
            logger.warning("overwriting registered tool %s", name)
        _TOOL_REGISTRY[name] = fn
        logger.debug("registered tool '%s'", name)
        return fn

    return deco


def get_tool(name: str) -> Callable[..., Any]:
    return _TOOL_REGISTRY[name]


def list_tools() -> list[str]:
    return sorted(_TOOL_REGISTRY)


@register_tool("lookup_gene")
def lookup_gene(name_or_id: str, bridge_path: Optional[str] = None) -> str:
    """Look up a gene name or dd_Smed_v6 / v4 id in the bridge table.

    Parameters
    ----------
    name_or_id
        Either a plain gene name or a dd_Smed_v6/dd_Smed_v4 id string.
    bridge_path
        Optional path to a bridge CSV. If omitted, the tool returns the
        query itself so an LLM can still reason about it. The NeuralTF
        project (Layer 9) injects its real bridge path when calling.
    """
    if not bridge_path:
        return json.dumps(
            {"ok": False, "note": "no bridge table available", "query": name_or_id}
        )
    p = Path(bridge_path)
    if not p.exists():
        return json.dumps(
            {"ok": False, "note": f"bridge not found at {bridge_path}", "query": name_or_id}
        )
    import pandas as pd

    df = pd.read_csv(p)
    hits: dict[str, Any] = {}
    if "gene_name" in df.columns:
        row = df.loc[df["gene_name"] == name_or_id]
    else:
        row = df.iloc[0:0]
    if row.empty and "v6_id" in df.columns:
        row = df.loc[df["v6_id"] == name_or_id]
    if row.empty and "v4_id" in df.columns:
        row = df.loc[df["v4_id"] == name_or_id]
    if row.empty:
        return json.dumps({"ok": False, "note": "no match", "query": name_or_id})
    rec = row.iloc[0].to_dict()
    return json.dumps({"ok": True, "gene": rec}, default=str)


@register_tool("summarize_candidates")
def summarize_candidates(rank_csv: str, top_n: int = 5) -> str:
    """Summarize the top-ranked candidate TFs from a CSV.

    The CSV must have at least columns ``gene_id``, ``tier``, and
    ``integrated_score`` (output produced by Layer 8B / Layer 9 ranking).
    Returns a JSON string with the top-N rows so an LLM assistant can
    write a one-paragraph "what to chase next" note.
    """
    p = Path(rank_csv)
    if not p.exists():
        return json.dumps({"ok": False, "note": f"{rank_csv} not found"})
    rows = []
    with open(p, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rows.append(r)
    if not rows:
        return json.dumps({"ok": False, "note": "no rows in CSV"})
    keep = rows[:top_n]
    return json.dumps({"ok": True, "top": keep, "n_total": len(rows)})
