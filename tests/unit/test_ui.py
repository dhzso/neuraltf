"""Smoke tests for bioforge.ui (Layer 10).

We test the pure-logic helpers under app.py directly so the tests don't
need a live Streamlit context.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from bioforge.ui.app import _load_workflow_yaml


def test_load_workflow_yaml_parses_steps(tmp_path: Path) -> None:
    yaml_text = """
description: Tiny workflow for testing.
steps:
  - id: greet
    target: upper
    inputs: {text: hi}
"""
    p = tmp_path / "wf.yaml"
    p.write_text(yaml_text)
    doc = _load_workflow_yaml(p)
    assert doc["description"].startswith("Tiny")
    assert doc["steps"][0]["id"] == "greet"


def test_app_main_module_imports_without_streamlit_run() -> None:
    # Plain import sanity check; streamlit is installed in the container so this
    # also exercises the real Streamlit import branch.
    import bioforge.ui.app as app
    assert callable(app.main)
    assert callable(app.render_run_page)
    assert callable(app.render_results_page)
    assert callable(app.render_assistant_page)
