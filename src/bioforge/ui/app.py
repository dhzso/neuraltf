"""Streamlit app — Run / Results / Assistant pages.

Designed to be importable without streamlit installed (for unit-testing
the pure logic that lives under the surface). When run via
``streamlit run src/bioforge/ui/app.py`` it builds the pages. The AI
assistant panel uses :func:`bioforge.ai.build_assistant`, which falls back
to :class:`StubAssistant` when no API key is configured.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

# Streamlit import is deferred so the package imports cleanly even if the
# streamlit extra isn't installed (for the UI's pure-logic helpers in tests).
def _st():
    import streamlit as st
    return st


def _load_workflow_yaml(path: Path) -> dict:
    import yaml
    return yaml.safe_load(path.read_text())


# ---------------------------------------------------------------------------
# Pages
# ---------------------------------------------------------------------------


def render_run_page() -> None:
    st = _st()
    st.subheader("Run a discovery workflow")
    workflow_dir = Path("projects/NeuralTF/workflows")
    if not workflow_dir.exists():
        st.warning("No project workflows found; create one via `bioforge projects create`.")
        return
    yamls = sorted(workflow_dir.glob("*.yaml")) + sorted(workflow_dir.glob("*.yml"))
    if not yamls:
        st.warning(f"No workflow YAMLs in {workflow_dir}")
        return
    selection = st.selectbox("Workflow", [p.name for p in yamls])
    datasets_input = st.text_area(
        "Datasets (one per line)",
        placeholder="Enter GEO/SRA accession, URL, or local file path...",
    )
    if st.button("Run workflow"):
        _run_workflow(workflow_dir / selection, datasets_input)


def _run_workflow(workflow_path: Path, datasets_input: str) -> None:
    import datetime as dt
    st = _st()
    from bioforge.workflow import WorkflowExecutor, WorkflowRun
    import bioforge.workflow.steps  # noqa: F401 — registers steps

    ts = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dir = Path("projects/NeuralTF/runs") / ts
    out_dir.mkdir(parents=True, exist_ok=True)
    progress = st.progress(0.0, text="Loading workflow...")
    log_box = st.empty()
    log_lines: list[str] = []

    def cb(step_id: str, target: str, duration: float) -> None:
        log_lines.append(f"  ✓ {step_id} ({target}) in {duration:.2f}s")
        log_box.markdown("\n".join(log_lines))

    run = WorkflowRun.from_yaml(workflow_path)
    executor = WorkflowExecutor(progress_cb=cb)
    try:
        outputs = executor.execute(run)
    except Exception as exc:  # noqa: BLE001
        st.error(f"Workflow failed: {exc}")
        return
    progress.progress(1.0, text="Done")
    (out_dir / "provenance.json").write_text(
        json.dumps(executor.provenance, indent=2), encoding="utf-8"
    )
    summary = {
        "workflow_yaml": str(workflow_path),
        "out_dir": str(out_dir),
        "n_steps": len(run.steps),
        "step_ids": [s.id for s in run.steps],
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    st.success(f"Workflow complete. Artifacts in {out_dir}")
    st.session_state["last_run_dir"] = str(out_dir)


def render_results_page() -> None:
    st = _st()
    st.subheader("Results")
    last = st.session_state.get("last_run_dir")
    runs_root = Path("projects/NeuralTF/runs")
    candidates = sorted([p for p in runs_root.glob("*") if p.is_dir()], reverse=True) if runs_root.exists() else []
    if last:
        candidates = [Path(last)] + [p for p in candidates if str(p) != last]
    if not candidates:
        st.info("No runs available yet. Click 'Run workflow' in the Run tab.")
        return
    run_path = st.selectbox("Run", [str(p) for p in candidates])
    run = Path(run_path)
    cards_md = run / "evidence_cards.md"
    rank_csv = run / "rank.csv"
    if rank_csv.exists():
        try:
            import pandas as pd
            df = pd.read_csv(rank_csv)
            st.dataframe(df)
        except Exception as exc:  # noqa: BLE001
            st.warning(f"Could not display rank CSV: {exc}")
    if cards_md.exists():
        st.markdown(cards_md.read_text(encoding="utf-8"))
    else:
        st.info(
            "No evidence cards were written in this run. The workflow YAML "
            "needs a `report.write_cards_md` step to surface findings here."
        )


def render_assistant_page() -> None:
    st = _st()
    st.subheader("AI Assistant")
    st.caption(
        "Configure by exporting `BIOFORGE_AI_API_KEY` (and optionally "
        "`BIOFORGE_AI_BASE_URL` / `BIOFORGE_AI_MODEL`) before launching Streamlit."
    )
    from bioforge.ai import build_assistant
    from bioforge.ai.assistant import ChatMessage
    assistant = build_assistant()
    st.info(f"Active assistant: **{assistant.name}**")
    if "chat_history" not in st.session_state:
        st.session_state["chat_history"] = []
    user_msg = st.chat_input("Ask about candidates, datasets, or the workflow...")
    if user_msg:
        st.session_state["chat_history"].append(("user", user_msg))
        try:
            resp = assistant.complete([ChatMessage(role="user", content=user_msg)])
        except Exception as exc:  # noqa: BLE001
            resp = type("R", (), {"content": f"[error] {exc}", "model": "error", "elapsed_s": 0.0})()
        st.session_state["chat_history"].append(("assistant", resp.content))
    for role, content in st.session_state["chat_history"]:
        with st.chat_message(role):
            st.markdown(content)


def main() -> None:
    """Streamlit entry point — launching with `streamlit run app.py` calls this."""
    st = _st()
    st.set_page_config(page_title="BioForge", layout="wide")
    st.title("BioForge")
    st.caption("AI-native planarian TF discovery workstation (NeuralTF thesis).")
    page = st.sidebar.selectbox("Page", ["Run", "Results", "Assistant"])
    if page == "Run":
        render_run_page()
    elif page == "Results":
        render_results_page()
    else:
        render_assistant_page()


if __name__ == "__main__":  # pragma: no cover — exercised via `streamlit run`
    if "--help" in sys.argv:
        print("Streamlit app. Launch with: streamlit run src/bioforge/ui/app.py")
        sys.exit(0)
    main()
