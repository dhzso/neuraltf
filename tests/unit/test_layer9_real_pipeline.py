"""Tests for `bioforge run` CLI command + new step registrations + loop keyword."""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from bioforge.cli.main import cli
from bioforge.evidence.schema import EvidenceSource
from bioforge.workflow import WorkflowExecutor, WorkflowRun, WorkflowStep
from bioforge.workflow.registry import StepRegistry as _SR


# ---------------------------------------------------------------------------
# loop: keyword in WorkflowStep / WorkflowExecutor
# ---------------------------------------------------------------------------


def _reg_with_loop_steps():
    reg = _SR()
    reg.register("echo", lambda x: {"value": x})
    reg.register("append_one", lambda x: x + 1)
    reg.register("triple", lambda x: {"result": [x, x, x]})
    return reg


def test_loop_runs_step_once_per_list_item() -> None:
    reg = _reg_with_loop_steps()
    ex = WorkflowExecutor(registry=reg)
    run = WorkflowRun(
        steps=[
            WorkflowStep(id="triple", target="triple", inputs={"x": 1}),
            WorkflowStep(id="echo_each", target="echo",
                         inputs={"x": "$item"}, loop="$triple.result"),
        ],
    )
    out = ex.execute(run)
    assert out["echo_each"]["results"][0]["value"] == 1
    assert out["echo_each"]["results"][1]["value"] == 1
    assert out["echo_each"]["results"][2]["value"] == 1
    # provenance records 1 (triple) + 3 (loop iterations)
    loop_provenance = [p for p in ex.provenance if p["step_id"] == "echo_each"]
    assert len(loop_provenance) == 3
    for p in loop_provenance:
        assert "stdout_tail" in p


def test_loop_with_literal_list_param() -> None:
    reg = _reg_with_loop_steps()
    ex = WorkflowExecutor(registry=reg)
    run = WorkflowRun(steps=[
        WorkflowStep(id="echo_each", target="echo",
                     inputs={"x": "$item"}, loop=["a", "b"]),
    ])
    out = ex.execute(run)
    assert [r["value"] for r in out["echo_each"]["results"]] == ["a", "b"]


def test_workflow_inputs_injected_via_extra_inputs() -> None:
    reg = _SR()
    reg.register("echo", lambda x: {"value": x})
    ex = WorkflowExecutor(registry=reg)
    run = WorkflowRun(steps=[
        WorkflowStep(id="ext", target="echo", inputs={"x": "$inputs.user_value"}),
    ])
    out = ex.execute(run, extra_inputs={"user_value": 42})
    assert out["ext"]["value"] == 42


# ---------------------------------------------------------------------------
# New workflow steps (load_bridge / score_per_atlas / combine_scores / ...)
# ---------------------------------------------------------------------------


def test_step_load_bridge_csv_round_trip(tmp_path: Path) -> None:
    import bioforge.workflow.steps as _steps  # noqa: F401
    df = pd.DataFrame({
        "gene_name": ["soxB"], "v6_id": ["v6_1"], "v4_id": ["v4_1"],
    })
    p = tmp_path / "bridge.csv"
    df.to_csv(p, index=False)

    from bioforge.workflow.registry import StepRegistry
    fn = StepRegistry.instance().get("evidence.load_bridge_csv")
    out = fn(path=str(p))
    assert out["bridge"].n_bridged == 1


def test_step_combine_scores_aggregates_by_gene() -> None:
    import bioforge.workflow.steps as _steps  # noqa: F401
    from bioforge.evidence.schema import EvidenceRecord, EvidenceSource
    from bioforge.workflow.registry import StepRegistry

    a_rec = EvidenceRecord(gene_id="g1", gene_name="soxB")
    a_rec.add_score(EvidenceSource.EXPRESSION, 0.8, note="atlas=f")
    b_rec = EvidenceRecord(gene_id="g1", gene_name="soxB")
    b_rec.add_score(EvidenceSource.EXPRESSION, 1.0, note="atlas=k")
    fn = StepRegistry.instance().get("evidence.combine_scores")
    out = fn(records_per_atlas=[
        {"records": [a_rec], "atlas_name": "fincher"},
        {"records": [b_rec], "atlas_name": "king"},
    ])
    assert len(out["records"]) == 1
    merged = out["records"][0]
    # max-of-expression = 1.0
    assert merged.scores[EvidenceSource.EXPRESSION] == 1.0
    # reproducibility = 2 atlases / 3 max = 2/3
    assert abs(merged.scores[EvidenceSource.REPRODUCIBILITY] - 2.0 / 3.0) < 1e-6
    assert out["atlas_membership"]["g1"] == ["fincher", "king"]


def test_step_add_rnai_stream_sets_score_when_in_table() -> None:
    import bioforge.workflow.steps as _steps  # noqa: F401
    from bioforge.evidence.schema import EvidenceRecord, EvidenceSource
    from bioforge.workflow.registry import StepRegistry

    rnai = pd.DataFrame({"fstf_rnai": ["soxB"], "marker": ["m1"]})
    recs = [EvidenceRecord(gene_id="g1", gene_name="soxB")]
    fn = StepRegistry.instance().get("evidence.add_rnai_stream")
    out = fn(records=recs, rnai_table=rnai)
    assert out["records"][0].scores[EvidenceSource.RNai] == 1.0


def test_step_add_rnai_stream_zero_when_not_in_table() -> None:
    import bioforge.workflow.steps as _steps  # noqa: F401
    from bioforge.evidence.schema import EvidenceRecord
    from bioforge.workflow.registry import StepRegistry

    rnai = pd.DataFrame({"fstf_rnai": ["soxB"], "marker": ["m1"]})
    recs = [EvidenceRecord(gene_id="g2", gene_name="novelTf")]
    fn = StepRegistry.instance().get("evidence.add_rnai_stream")
    out = fn(records=recs, rnai_table=rnai)
    assert out["records"][0].scores.get(EvidenceSource.RNai, 0.0) == 0.0


# ---------------------------------------------------------------------------
# CLI run --input flag + artifacts dir + ai_summary.md
# ---------------------------------------------------------------------------


def test_cli_run_help_includes_inputs_option() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "--help"])
    assert result.exit_code == 0
    assert "--input" in result.output


def test_cli_run_writes_artifacts_dir_and_provenance(tmp_path: Path) -> None:
    runner = CliRunner()
    wf = tmp_path / "wf.yaml"
    wf.write_text("""
description: Demo
steps:
  - id: r
    target: evidence.demo_rank
    inputs: {}
    params: {}
""")
    out_dir = tmp_path / "out"
    result = runner.invoke(cli, ["run", str(wf), "--out", str(out_dir)])
    assert result.exit_code == 0, result.output
    assert (out_dir / "provenance.json").exists()
    assert (out_dir / "summary.json").exists()
    assert (out_dir / "artifacts").is_dir()


def test_cli_run_with_input_kv_overrides_workflow_inputs(tmp_path: Path) -> None:
    reg = _SR()
    reg.register("echo_str", lambda x: {"result": x})
    # Use a temporary registry by directly invoking the executor (CLI is
    # always wired to the singleton registry; here we exercise the
    # extra_inputs pathway via the executor API).
    ex = WorkflowExecutor(registry=reg)
    run = WorkflowRun(steps=[WorkflowStep(id="e", target="echo_str",
                                          inputs={"x": "$inputs.disease"})])
    out = ex.execute(run, extra_inputs={"disease": "neural"})
    assert out["e"]["result"] == "neural"


# ---------------------------------------------------------------------------
# AI missing tool — inspect_anndata
# ---------------------------------------------------------------------------


def test_inspect_anndata_reads_summary(tmp_path: Path) -> None:
    import anndata as _ad
    _ad.settings.allow_write_nullable_strings = True
    _pd_set_string_python()
    obs = pd.DataFrame(index=np.array([f"c{i}" for i in range(3)], dtype=object))
    var = pd.DataFrame(index=np.array([f"g{i}" for i in range(3)], dtype=object))
    adata = ad.AnnData(X=np.eye(3, dtype=np.float32), obs=obs, var=var)
    adata.obs["leiden"] = ["0", "1", "0"]
    p = tmp_path / "x.h5ad"
    adata.write_h5ad(p)
    from bioforge.ai import inspect_anndata
    import json as _json
    payload = _json.loads(inspect_anndata(str(p)))
    assert payload["ok"] is True
    assert payload["n_obs"] == 3 and payload["n_vars"] == 3
    assert "leiden" in payload["obs_columns"]


def test_inspect_anndata_returns_false_for_missing_file(tmp_path: Path) -> None:
    import json as _json
    from bioforge.ai import inspect_anndata
    payload = _json.loads(inspect_anndata(str(tmp_path / "nope.h5ad")))
    assert payload["ok"] is False


def _pd_set_string_python():
    import pandas as _pd
    _pd.set_option("mode.string_storage", "python")


# ---------------------------------------------------------------------------
# WorkflowRun.from_yaml parses loop: field
# ---------------------------------------------------------------------------


def test_workflow_run_from_yaml_reads_loop_field(tmp_path: Path) -> None:
    wf = tmp_path / "wf.yaml"
    wf.write_text("""
steps:
  - id: e
    target: echo
    inputs: {x: $item}
    loop: $prev.list
""")
    run = WorkflowRun.from_yaml(wf)
    assert run.steps[0].loop == "$prev.list"


# ---------------------------------------------------------------------------
# Provenance captures stdout_tail
# ---------------------------------------------------------------------------


def test_provenance_captures_stdout_tail() -> None:
    reg = _SR()
    reg.register("printy", lambda: print("hello-from-step"))
    ex = WorkflowExecutor(registry=reg)
    run = WorkflowRun(steps=[WorkflowStep(id="p", target="printy")])
    ex.execute(run)
    assert "hello-from-step" in ex.provenance[0]["stdout_tail"]


# ---------------------------------------------------------------------------
# NeuralTF real.yaml parses cleanly (doesn't have to run; just structurally valid)
# ---------------------------------------------------------------------------


def test_real_pipeline_yaml_parses() -> None:
    wf = Path("projects/NeuralTF/workflows/real.yaml")
    run = WorkflowRun.from_yaml(wf)
    assert len(run.steps) >= 10
    ids = [s.id for s in run.steps]
    for required in ["ingest_fincher", "ingest_plass", "ingest_king",
                     "load_bridge", "score_per_atlas", "combine_scores",
                     "assign_tiers", "build_cards", "write_cards",
                     "summarise"]:
        assert required in ids, f"missing required step: {required}"


def test_discover_neural_tfs_yaml_filename_is_present() -> None:
    # README documents `bioforge run projects/NeuralTF/workflows/discover_neural_tfs.yaml`
    assert Path("projects/NeuralTF/workflows/discover_neural_tfs.yaml").exists()
