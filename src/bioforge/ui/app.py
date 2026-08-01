"""Streamlit app for the NeuralTF pipeline.

Three pages:
  - **Run**     : one-click execution of `NeuralTFPipeline` with live progress
                  and parameters (subsample size, output directory).
  - **Results** : browse `rank.csv` / `rank_neural.csv` from any run, filter
                  by tier, search by gene, and read the markdown evidence
                  cards.
  - **Assistant**: chat with the AI assistant (StubAssistant when no API
                  key is configured) for follow-up queries about candidates.

The app is importable without streamlit installed (the `_st()` indirection
defers the import so the pure-logic helpers under the surface can be unit
tested). Launch with `bioforge ui` or `streamlit run src/bioforge/ui/app.py`.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

# Streamlit import is deferred so the package imports cleanly even if the
# streamlit extra isn't installed (the pure-logic helpers stay unit-testable).
def _st():
    import streamlit as st
    return st


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _repo_root() -> Path:
    """Return the repo root: BIOFORGE_REPO_ROOT env (set by `bioforge ui`)
    if present, else the cwd."""
    import os
    return Path(os.environ.get("BIOFORGE_REPO_ROOT", Path.cwd()))


def _runs_root() -> Path:
    return _repo_root() / "projects" / "NeuralTF" / "runs"


def _list_runs() -> list[Path]:
    root = _runs_root()
    if not root.exists():
        return []
    return sorted(
        [p for p in root.iterdir() if p.is_dir()],
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )


def _pipeline_status() -> dict:
    """Return a short status snapshot about data files and prior runs."""
    root = _repo_root()
    proc = root / "datasets" / "processed"
    fincher = proc / "fincher_subsample.h5ad"
    plass = proc / "plass_v6.h5ad"
    king_dir = root / "datasets" / "raw" / "Supplementary_Data_ King_2024"
    return {
        "fincher_h5ad": fincher.exists(),
        "fincher_path": str(fincher),
        "plass_h5ad": plass.exists(),
        "plass_path": str(plass),
        "bridge_csv": (root / "projects" / "NeuralTF" / "data" / "bridge.csv").exists(),
        "king_atlas_tsv": (root / "projects" / "NeuralTF" / "data" / "king_atlas.tsv").exists(),
        "king_dir_exists": king_dir.exists(),
        "n_runs": len(_list_runs()),
    }


# ---------------------------------------------------------------------------
# Pages
# ---------------------------------------------------------------------------


def render_run_page() -> None:
    st = _st()
    st.subheader("Run the NeuralTF pipeline")

    status = _pipeline_status()
    ready = status["fincher_h5ad"] and status["plass_h5ad"]
    if not ready:
        st.warning("Required processed h5ad files are missing.")
        st.markdown(
            f"- `fincher_subsample.h5ad`: {'OK' if status['fincher_h5ad'] else 'MISSING'} "
            f"({status['fincher_path']})\n"
            f"- `plass_v6.h5ad`: {'OK' if status['plass_h5ad'] else 'MISSING'} "
            f"({status['plass_path']})\n\n"
            "Build them from raw GEO downloads:\n"
            "```\n"
            "python scripts/convert_fincher.py\n"
            "python scripts/consolidate_plass.py\n"
            "```"
        )
        st.caption(
            "Sources: GEO [GSE111764](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111764) "
            "(Fincher 2018) and [GSE109226](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109226) "
            "(Plass 2018)."
        )
        if not status["king_dir_exists"]:
            st.error(
                "King 2024 supplementary xlsx (mmc4-mmc7) is also missing. Place under "
                "`datasets/raw/Supplementary_Data_ King_2024/`."
            )
        return

    st.success("All inputs present. Ready to run.")

    with st.form("pipeline_params"):
        col1, col2 = st.columns(2)
        with col1:
            subsample = st.number_input(
                "Cells per atlas (subsample)",
                min_value=0,
                max_value=40000,
                value=10000,
                step=1000,
                help="Set to 0 for no subsampling (use all available cells).",
            )
        with col2:
            out_dir = st.text_input(
                "Output directory (relative to repo root)",
                value="projects/NeuralTF/runs/pipeline_run",
                help="Where rank.csv, rank_neural.csv and evidence_cards.md are written.",
            )
        submitted = st.form_submit_button("Run pipeline", type="primary")

    if not submitted:
        # Show existing runs as a hint
        runs = _list_runs()
        if runs:
            st.markdown(f"**{len(runs)} existing run(s)** — see the **Results** page.")
        return

    from bioforge.projects.neuraltf.pipeline import NeuralTFPipeline

    root = _repo_root()
    out_path = root / out_dir
    st.info(f"Starting pipeline. Output: `{out_path}`")
    progress = st.progress(0.0, text="Initializing...")
    log_box = st.empty()
    log_lines: list[str] = []

    # Capture print() output of the pipeline by redirecting to our log box.
    import io
    import contextlib

    class _StreamlitLog(io.StringIO):
        def write(self, s: str) -> int:
            log_lines.append(s.rstrip())
            log_box.markdown(
                "```\n" + "\n".join(log_lines[-200:]) + "\n```"
            )
            return len(s)

    # Build pipeline. Override default paths to use the user's repo root.
    pipe = NeuralTFPipeline(
        data_root=root,
        out_dir=out_path,
        subsample=int(subsample) if subsample > 0 else None,
    )

    # Step atomically through the pipeline so we can update the progress bar.
    steps = [
        ("load_datasets", "1/8 Loading datasets"),
        ("load_reference_tables", "2/8 Reference tables"),
        ("run_qc", "4/8 QC + clustering"),
        ("score_atlases", "5/8 Atlas DE scoring"),
        ("integrate_king_atlas", "6/8 King neural seed"),
        ("integrate_rnai", "7/8 RNAi evidence"),
        ("integrate_correlations", "7.5/8 Correlations"),
        ("assign_reproducibility", "8/8 Reproducibility"),
        ("write_outputs", "Writing outputs"),
    ]
    log_buf = _StreamlitLog()
    try:
        with contextlib.redirect_stdout(log_buf):
            for i, (method_name, label) in enumerate(steps, start=1):
                progress.progress(i / len(steps), text=label)
                getattr(pipe, method_name)()
        progress.progress(1.0, text="Done")
        st.success(f"Pipeline complete. Artifacts in `{out_path}`")
        st.session_state["last_run_dir"] = str(out_path)
        # Provide quick links
        for fname in ("rank.csv", "rank_neural.csv", "evidence_cards.md", "pipeline_results.json"):
            p = out_path / fname
            if p.exists():
                st.markdown(f"- `{p.relative_to(root)}` ({p.stat().st_size :,} bytes)")
        st.button("View results →", on_click=lambda: st.session_state.update(page="Results"))
    except Exception as exc:  # noqa: BLE001
        progress.progress(0.0, text="Failed")
        st.error(f"Pipeline failed: {exc}")
        with st.expander("Log"):
            st.code("\n".join(log_lines[-200:]))


def render_results_page() -> None:
    st = _st()
    st.subheader("Results")

    import pandas as pd

    runs = _list_runs()
    if not runs:
        st.info("No runs yet. Go to the **Run** page.")
        return

    selected = st.selectbox(
        "Run",
        [str(p.relative_to(_repo_root())) for p in runs],
        index=0,
    )
    run = _repo_root() / selected
    rank_csv = run / "rank.csv"
    neural_csv = run / "rank_neural.csv"
    cards_md = run / "evidence_cards.md"
    json_path = run / "pipeline_results.json"

    # Tabbed view of the two rankings + cards + run metadata
    tab_rank, tab_neural, tab_cards, tab_meta = st.tabs(
        ["All candidates", "Neural-filtered", "Evidence cards", "Run metadata"]
    )

    with tab_rank:
        if rank_csv.exists():
            df = pd.read_csv(rank_csv)
            st.caption(f"{len(df)} rows × {len(df.columns)} columns")
            _render_rank_table(df, key_prefix="all")
        else:
            st.warning("rank.csv not found in this run.")

    with tab_neural:
        if neural_csv.exists():
            df = pd.read_csv(neural_csv)
            st.caption(f"{len(df)} rows × {len(df.columns)} columns")
            _render_rank_table(df, key_prefix="neural")
            if "tier" in df.columns:
                st.bar_chart(df["tier"].value_counts())
        else:
            st.warning("rank_neural.csv not found in this run.")

    with tab_cards:
        if cards_md.exists():
            md = cards_md.read_text(encoding="utf-8")
            # Split on candidate delimiter (## headings) into a searchable list
            candidates = [c for c in md.split("\n## ") if c.strip()]
            if candidates and not candidates[0].startswith("#"):
                candidates = candidates[1:]  # drop preamble before first ##
            if candidates:
                labels = [c.split("\n", 1)[0].strip() for c in candidates]
                choice = st.selectbox("Candidate", labels)
                idx = labels.index(choice)
                st.markdown("## " + candidates[idx])
            else:
                st.markdown(md)
        else:
            st.warning("evidence_cards.md not found in this run.")

    with tab_meta:
        if json_path.exists():
            try:
                meta = json.loads(json_path.read_text(encoding="utf-8"))
                st.json(meta)
            except Exception as exc:  # noqa: BLE001
                st.warning(f"Could not parse {json_path.name}: {exc}")
        else:
            st.info("No pipeline_results.json in this run.")


def _render_rank_table(df, *, key_prefix: str = "rank") -> None:
    """Render a rank dataframe with a search + tier filter.

    ``key_prefix`` must be unique per caller so Streamlit's widget IDs don't
    collide when the function is invoked twice on the same page (e.g. one
    for the all-candidates tab and one for the neural-filtered tab).
    """
    st = _st()
    import pandas as pd

    # Search box for gene name / id
    search = st.text_input("Search gene", "", key=f"{key_prefix}_search")
    if search:
        mask = pd.Series(False, index=df.index)
        for col in ("gene_name", "gene_id", "gene_id_v6"):
            if col in df.columns:
                mask |= df[col].astype(str).str.contains(search, case=False, na=False)
        df = df[mask]

    # Tier filter
    if "tier" in df.columns:
        tiers = ["(all)"] + sorted(df["tier"].dropna().unique().tolist())
        choice = st.selectbox("Tier", tiers, key=f"{key_prefix}_tier")
        if choice != "(all)":
            df = df[df["tier"] == choice]

    # Show key columns upfront, full table under expander
    key_cols = [c for c in ("rank", "gene_name", "gene_id", "tier", "score", "proof_status") if c in df.columns]
    st.dataframe(df[key_cols] if key_cols else df, use_container_width=True)
    with st.expander("Show all columns"):
        st.dataframe(df, use_container_width=True)


def render_assistant_page() -> None:
    st = _st()
    st.subheader("AI Assistant")
    st.caption(
        "Configure by exporting `BIOFORGE_AI_API_KEY` (and optionally "
        "`BIOFORGE_AI_BASE_URL` / `BIOFORGE_AI_MODEL`) before launching."
    )

    from bioforge.ai import build_assistant
    from bioforge.ai.assistant import ChatMessage
    from bioforge.ai.tools import get_tool, list_tools

    assistant = build_assistant()
    st.info(
        f"Active assistant: **{assistant.name}** · "
        f"Available tools: {', '.join(list_tools()) or '(none)'}"
    )
    if "chat_history" not in st.session_state:
        st.session_state["chat_history"] = []

    user_msg = st.chat_input("Ask about candidates, datasets, or the workflow...")
    if user_msg:
        st.session_state["chat_history"].append(("user", user_msg))
        reply_content: str
        try:
            if user_msg.startswith("tool:"):
                tokens = user_msg[len("tool:"):].strip().split(maxsplit=1)
                if len(tokens) >= 2:
                    tool_name, args = tokens[0], tokens[1].strip()
                    if tool_name in list_tools():
                        fn = get_tool(tool_name)
                        result = fn(*args.split("|"))
                        reply_content = f"Tool `{tool_name}` returned:\n```\n{result}\n```"
                    else:
                        reply_content = f"Unknown tool `{tool_name}`. Available: {list_tools()}"
                else:
                    reply_content = "Provide tool args after `tool: name args`."
            else:
                resp = assistant.complete([ChatMessage(role="user", content=user_msg)])
                reply_content = resp.content
        except Exception as exc:  # noqa: BLE001
            reply_content = f"[error] {exc}"
        st.session_state["chat_history"].append(("assistant", reply_content))

    for role, content in st.session_state["chat_history"]:
        with st.chat_message(role):
            st.markdown(content)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Streamlit entry point — `streamlit run app.py` invokes this."""
    st = _st()
    st.set_page_config(
        page_title="BioForge · NeuralTF",
        page_icon="🧠",
        layout="wide",
    )
    st.title("BioForge · NeuralTF")
    st.caption(
        "Planarian neural-fate-specific transcription factor candidate "
        "discovery (3 atlases, 8 evidence streams)."
    )
    page = st.sidebar.selectbox("Page", ["Run", "Results", "Assistant"])
    if page == "Run":
        render_run_page()
    elif page == "Results":
        render_results_page()
    else:
        render_assistant_page()


if __name__ == "__main__":  # pragma: no cover — exercised via `streamlit run`
    if "--help" in sys.argv:
        print("Streamlit app. Launch with: bioforge ui  or  streamlit run src/bioforge/ui/app.py")
        sys.exit(0)
    main()
