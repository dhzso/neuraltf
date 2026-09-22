"""Workflow dataclasses and executor with provenance capture."""
from __future__ import annotations

import contextlib
import copy
import datetime as dt
import hashlib
import io
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
    # When `loop` is set, the step is repeated once per item in the value
    # the loop resolves to. The loop value can be: a list literal (in
    # params/inputs), a `$step.output` ref that resolves to a list, or the
    # special string `$inputs.<key>` to pull a list from WorkflowRun.inputs.
    loop: Optional[str | list] = None


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
                loop=s.get("loop"),
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

    def _resolve_ref(self, ref: Any, outputs: dict[str, dict[str, Any]],
                     workflow_inputs: dict[str, Any] | None = None) -> Any:
        """Resolve a `{"$step": ..., "$output": ...}` or "$step.output" ref.

        Also resolves `$inputs.<key>` references against the workflow-level
        inputs (so the CLI/UI can inject datasets per-run).
        """
        if isinstance(ref, dict) and "$step" in ref:
            step_id = ref["$step"]
            output_key = ref.get("$output") or _first_key_after(ref, "$step")
            return outputs.get(step_id, {}).get(output_key)
        if isinstance(ref, str) and ref.startswith("$"):
            ref_inner = ref[1:]
            if ref_inner.startswith("inputs."):
                key = ref_inner.split(".", 1)[1]
                return (workflow_inputs or {}).get(key)
            if "." in ref_inner:
                step_id, output_key = ref_inner.split(".", 1)
            else:
                step_id = ref_inner
                output_key = "result"
            return outputs.get(step_id, {}).get(output_key)
        return ref

    def _resolve_loop(self, step: WorkflowStep,
                      outputs: dict[str, dict[str, Any]],
                      workflow_inputs: dict[str, Any]) -> Optional[list]:
        """Resolve `step.loop` into a concrete list (or None)."""
        if step.loop is None:
            return None
        val = step.loop
        if isinstance(val, list):
            return val
        if isinstance(val, str):
            resolved = self._resolve_ref(val, outputs, workflow_inputs)
            if isinstance(resolved, list):
                return resolved
        return None

    def execute(self, run: WorkflowRun,
                extra_inputs: Optional[dict[str, Any]] = None) -> dict[str, dict[str, Any]]:
        '''Run all steps, returning a dict of step-id → outputs.

        `extra_inputs` are merged on top of WorkflowRun.inputs so the caller
        (CLI / UI) can override per-run values like datasets.
        '''
        outputs: dict[str, dict[str, Any]] = {}
        workflow_inputs = dict(run.inputs)
        if extra_inputs:
            workflow_inputs.update(extra_inputs)

        for step in run.steps:
            loop_items = self._resolve_loop(step, outputs, workflow_inputs)
            if loop_items is not None:
                step_results = []
                # If the per-iteration result is dict-valued, we collect the
                # list into outputs[step.id] under the "results" key; if
                # each iteration returns one AnnData-like object, the list
                # itself goes under "adatas" (or a loop_name picked from
                # items, when items are strings).
                for item in loop_items:
                    iter_inputs = {}
                    for k, v in step.inputs.items():
                        iter_inputs[k] = (
                            item if v == "$item" else
                            self._resolve_ref(v, outputs, workflow_inputs))
                    iter_params = copy.deepcopy(step.params)
                    for k, v in iter_params.items():
                        if isinstance(v, (dict, str)):
                            iter_params[k] = (
                                item if v == "$item" else
                                self._resolve_ref(v, outputs, workflow_inputs))
                    target = self.registry.get(step.target)
                    buf = io.StringIO()
                    t0 = dt.datetime.now(dt.timezone.utc)
                    try:
                        with contextlib.redirect_stdout(buf):
                            iter_out = target(**iter_inputs, **iter_params)
                    except Exception as exc:
                        logger.exception("step '%s' (loop item=%r) failed", step.id, item)
                        self.provenance.append({
                            "step_id": step.id, "target": step.target,
                            "started_at": t0.isoformat(), "duration_s": 0.0,
                            "error": str(exc), "input_hash": _hash_value(iter_inputs),
                            "loop_item": _hash_value(item),
                        })
                        raise
                    if not isinstance(iter_out, dict):
                        iter_out = {"result": iter_out}
                    step_results.append(iter_out)
                    duration = (dt.datetime.now(dt.timezone.utc) - t0).total_seconds()
                    self.provenance.append({
                        "step_id": step.id, "target": step.target,
                        "started_at": t0.isoformat(), "duration_s": duration,
                        "error": None, "input_hash": _hash_value(iter_inputs),
                        "outputs": list(iter_out.keys()),
                        "loop_item": _hash_value(item),
                        "stdout_tail": buf.getvalue()[-512:],
                    })
                outputs[step.id] = {"results": step_results}
                logger.info("step '%s' (%s) looped %d times", step.id, step.target, len(step_results))
                if self.progress_cb:
                    self.progress_cb(step.id, step.target, 0.0)
                continue

            t0 = dt.datetime.now(dt.timezone.utc)
            input_values: dict[str, Any] = {}
            params = copy.deepcopy(step.params)
            for k, v in step.inputs.items():
                input_values[k] = self._resolve_ref(v, outputs, workflow_inputs)
            for k, v in params.items():
                if isinstance(v, (dict, str)):
                    params[k] = self._resolve_ref(v, outputs, workflow_inputs)

            target = self.registry.get(step.target)
            buf = io.StringIO()
            try:
                with contextlib.redirect_stdout(buf):
                    fn_inputs = {**input_values}
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
                    "stdout_tail": buf.getvalue()[-512:],
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
                "stdout_tail": buf.getvalue()[-512:],
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
