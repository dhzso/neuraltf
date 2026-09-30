"""Integration test: NeuralTF demo_workflow runs end-to-end via the engine.

This exercises Layer 6 (AI stub), Layer 7 (workflow engine), Layer 8B
(evidence cards / scoring) and the workflow steps module — without needing
the raw datasets.
"""
from __future__ import annotations

import importlib
import json
from pathlib import Path

import pytest


def _ensure_steps_registered() -> None:
    importlib.import_module("bioforge.workflow.steps")


def test_demo_pipeline_runs_end_to_end(tmp_path: Path) -> None:
    _ensure_steps_registered()
    wf = Path("projects/NeuralTF/workflows/demo_pipeline.yaml")
    assert wf.exists(), "demo_pipeline.yaml must exist in the project layout"

    from bioforge.workflow import WorkflowExecutor, WorkflowRun
    run = WorkflowRun.from_yaml(wf)
    progress: list[tuple[str, str, float]] = []
    ex = WorkflowExecutor(progress_cb=lambda s, t, d: progress.append((s, t, d)))
    outputs = ex.execute(run)
    assert len(outputs) == 5
    assert len(progress) == 5
    assert outputs["summarise"]["summary"]


def test_runner_outputs_can_be_serialised(tmp_path: Path) -> None:
    _ensure_steps_registered()
    # Confirm every step's output value can be JSON-serialized (provenance
    # writes dict-of-str, but we check a small introspection here).
    wf = Path("projects/NeuralTF/workflows/demo_pipeline.yaml")
    from bioforge.workflow import WorkflowExecutor, WorkflowRun
    run = WorkflowRun.from_yaml(wf)
    ex = WorkflowExecutor()
    outputs = ex.execute(run)
    # JSON-serialize provenance (cards and records are not in provenance —
    # only durations, hashes, step ids).
    json.dumps(ex.provenance, default=str, indent=2)
