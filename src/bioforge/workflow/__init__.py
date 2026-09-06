"""BioForge Workflow Engine (Layer 7).

Declarative YAML workflows + a small Python-callables registry + a
dependency-resolving executor with provenance capture. Workflows are the
single end-to-end entry point used by both the CLI and the Streamlit UI.

The engine is deliberately tiny (no scheduler, no retries, no backpressure)
because we already have all the heavyweight compute in the omics layer; the
workflow engine's job is just to sequence steps, bind values, and record
what ran so a thesis run is reproducible.
"""
from bioforge.workflow.engine import WorkflowExecutor, WorkflowRun, WorkflowStep
from bioforge.workflow.registry import StepRegistry, register

__all__ = [
    "StepRegistry",
    "register",
    "WorkflowExecutor",
    "WorkflowRun",
    "WorkflowStep",
]
