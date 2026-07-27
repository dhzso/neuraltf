"""Workflow dataclasses and executor with provenance capture."""
from __future__ import annotations

import copy
import datetime as dt
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import yaml

from bioforge.core.logging import get_logger
from bioforge.workflow.registry import StepRegistry

logger = get_logger("workflow.engine")


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class WorkflowStep:
    id: str
    target: str
    inputs: dict[str, Any] = field(default_factory=dict)
    params: dict[str, Any] = field(default_factory=dict)
    description: str = ""


@dataclass
class StepResult:
    step_id: str
    outputs: dict[str, Any]
    duration_s: float
    error: Optional[str] = None


@dataclass
class WorkflowRun:
    steps: list[WorkflowStep]
    inputs: dict[str, Any] = field(default_factory=dict)
    description: str = ""

    @classmethod
    def from_yaml(cls, path: str | Path) -> "WorkflowRun":
        with open(path) as fh:
            doc = yaml.safe_load(fh)
        steps = [
            WorkflowStep(
                id=s["id"],
                target=s["target"],
                inputs=s.get("inputs", {}),
                params=s.get("params", {}),
                description=s.get("description", ""),
            )
            for s in doc.get("steps", [])
        ]
        return cls(
            steps=steps,
            inputs=doc.get("inputs", {}),
            description=doc.get("description", ""),
        )


# ---------------------------------------------------------------------------
# Executor
# ---------------------------------------------------------------------------


class WorkflowExecutor:
    """Runs a :class:`WorkflowRun`, binding values step by step.

    The executor exposes a tiny public API: :meth:`execute` returns a dict
    of step-id → outputs. The dict is stored on ``context`` so the caller
    (CLI / UI) can inspect provenance afterwards.

    References in YAML take two forms:
    - ``{"$step": "subcluster_qc", "$output": "adata"}`` (escaped as refs)
    - a Python dict string-of-shape `step_id.output_key` is also accepted
      when it's literally a plain string scalar starting with ``"$"``.
    """

    def __init__(
        self,
        registry: StepRegistry | None = None,
        progress_cb: Optional[Any] = None,
    ) -> None:
        self.registry = registry or StepRegistry.instance()
        self.progress_cb = progress_cb
        self.provenance: list[dict[str, Any]] = []

    def _resolve_ref(self, ref: Any, outputs: dict[str, dict[str, Any]]) -> Any:
        """Resolve a `{"$step": ..., "$output": ...}` or "$step.output" ref."""
        if isinstance(ref, dict) and "$step" in ref:
            step_id = ref["$step"]
            output_key = ref.get("$output") or _first_key_after(ref, "$step")
            return outputs.get(step_id, {}).get(output_key)
        if isinstance(ref, str) and ref.startswith("$"):
            ref_inner = ref[1:]
            if "." in ref_inner:
                step_id, output_key = ref_inner.split(".", 1)
            else:
                step_id = ref_inner
                output_key = "result"
            return outputs.get(step_id, {}).get(output_key)
        return ref

    def execute(self, run: WorkflowRun) -> dict[str, dict[str, Any]]:
        outputs: dict[str, dict[str, Any]] = {}
        workflow_inputs = dict(run.inputs)

        for step in run.steps:
            t0 = dt.datetime.now(dt.timezone.utc)
            input_values: dict[str, Any] = {}
            params = copy.deepcopy(step.params)
            for k, v in step.inputs.items():
                input_values[k] = self._resolve_ref(v, outputs)
            for k, v in params.items():
                if isinstance(v, (dict, str)):
                    params[k] = self._resolve_ref(v, outputs)

            target = self.registry.get(step.target)
            try:
                fn_inputs = {**input_values}
                # Workflow-level inputs become kwargs if they are not already
                # supplied via step.inputs ( .= stratifies explicit over env).
                step_outputs = target(**fn_inputs, **params)
            except Exception as exc:
                logger.exception("step '%s' failed", step.id)
                err = str(exc)
                self.provenance.append({
                    "step_id": step.id,
                    "target": step.target,
                    "started_at": t0.isoformat(),
                    "duration_s": 0.0,
                    "error": err,
                    "input_hash": _hash_value(input_values),
                })
                raise

            if not isinstance(step_outputs, dict):
                step_outputs = {"result": step_outputs}

            outputs[step.id] = step_outputs
            duration = (dt.datetime.now(dt.timezone.utc) - t0).total_seconds()
            self.provenance.append({
                "step_id": step.id,
                "target": step.target,
                "started_at": t0.isoformat(),
                "duration_s": duration,
                "error": None,
                "input_hash": _hash_value(input_values),
                "outputs": list(step_outputs.keys()),
            })
            logger.info("step '%s' (%s) done in %.3fs", step.id, step.target, duration)
            if self.progress_cb:
                self.progress_cb(step.id, step.target, duration)
        return outputs


def _first_key_after(d: dict, after: str) -> str:
    """Return the first key that isn't `after` from dict `d`; fallback 'result'."""
    for k in d:
        if k != after:
            return k
    return "result"


def _hash_value(value: Any) -> str:
    try:
        return hashlib.sha1(json.dumps(value, default=str).encode()).hexdigest()[:12]
    except Exception:
        return "unhashable"
