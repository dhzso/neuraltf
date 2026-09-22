"""Unit tests for bioforge.workflow (Layer 7)."""
from __future__ import annotations

from pathlib import Path

import pytest

from bioforge.workflow import (
    StepRegistry,
    WorkflowExecutor,
    WorkflowRun,
    WorkflowStep,
    register,
)
from bioforge.workflow.registry import StepRegistry as _SR


def _fresh_registry() -> StepRegistry:
    # Each test uses a brand new registry so global registrations don't leak.
    return _SR()


# ---------------------------------------------------------------------------
# StepRegistry
# ---------------------------------------------------------------------------


def test_step_registry_register_and_get() -> None:
    reg = _fresh_registry()
    reg.register("hello", lambda name: f"hi {name}")
    assert reg.get("hello")("deepa") == "hi deepa"
    assert "hello" in reg.known()


def test_step_registry_unknown_raises_keyerror() -> None:
    reg = _fresh_registry()
    with pytest.raises(KeyError):
        reg.get("does_not_exist")


def test_step_registry_schema_round_trip() -> None:
    reg = _fresh_registry()
    reg.register("s", lambda x: x, inputs=["x"], outputs=["y"], description="d")
    assert reg.schema("s")["inputs"] == ["x"]
    assert reg.schema("s")["outputs"] == ["y"]


# ---------------------------------------------------------------------------
# WorkflowExecutor
# ---------------------------------------------------------------------------


def test_executor_runs_simple_linear_workflow() -> None:
    reg = _fresh_registry()
    reg.register("upper", lambda text: text.upper())
    reg.register("exclaim", lambda text: f"{text}!")
    ex = WorkflowExecutor(registry=reg)
    run = WorkflowRun(
        steps=[
            WorkflowStep(id="up", target="upper", inputs={"text": "hello"}),
            WorkflowStep(id="shout", target="exclaim", inputs={"text": {"$step": "up", "$output": "result"}}),
        ],
    )
    out = ex.execute(run)
    assert out["up"]["result"] == "HELLO"
    assert out["shout"]["result"] == "HELLO!"
    assert len(ex.provenance) == 2
    assert all("duration_s" in p and "input_hash" in p for p in ex.provenance)


def test_executor_supports_dollar_string_refs() -> None:
    reg = _fresh_registry()
    reg.register("add_one", lambda x: x + 1)
    reg.register("double", lambda x: x * 2)
    ex = WorkflowExecutor(registry=reg)
    run = WorkflowRun(
        steps=[
            WorkflowStep(id="base", target="add_one", inputs={"x": 1}),
            WorkflowStep(id="dbl", target="double", inputs={"x": "$base.result"}),
        ],
    )
    out = ex.execute(run)
    assert out["base"]["result"] == 2
    assert out["dbl"]["result"] == 4


def test_executor_propagates_errors_and_records_provenance() -> None:
    reg = _fresh_registry()

    def boom(x):
        raise ValueError("kaboom")
    reg.register("boom", boom)
    ex = WorkflowExecutor(registry=reg)
    run = WorkflowRun(steps=[WorkflowStep(id="b", target="boom", inputs={"x": 1})])
    with pytest.raises(ValueError):
        ex.execute(run)
    assert len(ex.provenance) == 1
    assert "kaboom" in ex.provenance[0]["error"]


def test_executor_progress_callback_invoked() -> None:
    reg = _fresh_registry()
    reg.register("passthru", lambda x: x)
    seen = []
    ex = WorkflowExecutor(registry=reg, progress_cb=lambda s, t, d: seen.append(s))
    run = WorkflowRun(steps=[WorkflowStep(id="p", target="passthru", inputs={"x": 1})])
    ex.execute(run)
    assert seen == ["p"]


def test_workflow_run_from_yaml(tmp_path: Path) -> None:
    yaml_text = """
description: A tiny test workflow.
inputs: {}
steps:
  - id: greet
    target: upper
    inputs:
      text: hi
    params: {}
"""
    p = tmp_path / "wf.yaml"
    p.write_text(yaml_text)
    run = WorkflowRun.from_yaml(p)
    assert run.description.startswith("A tiny")
    assert run.steps[0].id == "greet"
    assert run.steps[0].target == "upper"
    assert run.steps[0].inputs == {"text": "hi"}


# ---------------------------------------------------------------------------
# Built-in registry default usage (smoke — the global registry singleton)
# ---------------------------------------------------------------------------


def test_register_decorator_attaches_to_singleton() -> None:
    @register("test_decorator_temp")
    def fn():
        return 42
    from bioforge.workflow.registry import StepRegistry
    assert StepRegistry.instance().get("test_decorator_temp")() == 42
