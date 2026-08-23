"""Streamlit app for the NeuralTF pipeline.

Four pages:
  - **Run**          : one-click execution of `NeuralTFPipeline` with live
                       progress and parameters (subsample size, output dir).
  - **Results**      : browse `rank.csv` / `rank_neural.csv` from any run,
                  filter by tier, search by gene, and read the markdown
                  evidence cards.
  - **Prioritization** : dual-track shortlist (`top10_neural_tfs_prioritized.csv`)
                  + `candidate_summary_report.md`.
  - **Assistant**   : chat with the AI assistant (StubAssistant when no API
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


def plt_close(fig) -> None:
    """Close a matplotlib Figure to free its memory after Streamlit renders
    it. Safe to call even if matplotlib isn't installed at import time."""
    try:
        import matplotlib.pyplot as plt
        plt.close(fig)
    except Exception:  # noqa: BLE001 - best-effort cleanup
        pass


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
    king_mmc = {
        f"mmc{i}": any(king_dir.glob(f"*mmc{i}*.xlsx"))
        if king_dir.exists() else False
        for i in (4, 5, 6, 7)
    }
    return {
        "fincher_h5ad": fincher.exists(),
        "fincher_path": str(fincher),
        "plass_h5ad": plass.exists(),
        "plass_path": str(plass),
        "bridge_csv": (root / "projects" / "NeuralTF" / "data" / "bridge.csv").exists(),
        "king_atlas_tsv": (root / "projects" / "NeuralTF" / "data" / "king_atlas.tsv").exists(),
        "king_dir_exists": king_dir.exists(),
        "king_mmc4_xlsx": king_mmc["mmc4"],
        "king_mmc5_xlsx": king_mmc["mmc5"],
        "king_mmc6_xlsx": king_mmc["mmc6"],
        "king_xlsx_ready": king_mmc["mmc4"] and king_mmc["mmc5"] and king_mmc["mmc6"],
        "n_runs": len(_list_runs()),
    }


def _downstream_readiness(root: Path) -> dict:
    """Which post-pipeline generation steps have their inputs available.

    ``figures`` needs only the completed run's rank.csv; ``prioritization``
    needs the King 2024 xlsx tables + the PlanMine parquet; ``go_figures``
    needs go.obo plus the shortlist CSV produced by prioritization.
    """
    king_dir = root / "datasets" / "raw" / "Supplementary_Data_ King_2024"
    mmc4 = any(king_dir.glob("*mmc4*.xlsx")) if king_dir.exists() else False
    mmc5 = any(king_dir.glob("*mmc5*.xlsx")) if king_dir.exists() else False
    return {
        "figures": True,
        "prioritization": (
            (root / "datasets" / "processed" / "planmine_annotations.parquet").exists()
            and mmc4
            and mmc5
        ),
        "go_figures": (root / "datasets" / "raw" / "go.obo").exists(),
    }


def _run_step(root: Path, st, log_lines: list[str], log_box,
              argv: list[str], step: str, ready: bool, why: str) -> bool:
    """Run one CLI step as a subprocess, teeing output into the log box.

    Returns True on success; skipped/missing-input and failures return False
    with a caption (never raises). Used by both the post-run pipeline steps
    and the data-rebuild steps below.
    """
    import subprocess
    import sys as _sys

    def _emit() -> None:
        log_box.markdown("```\n" + "\n".join(log_lines[-200:]) + "\n```")

    if not ready:
        st.caption(f"`{step}` skipped — {why}")
        return False
    log_lines.append(f"[run] python {' '.join(argv)}")
    _emit()
    try:
        res = subprocess.run(
            [_sys.executable, *argv], cwd=root, capture_output=True,
            text=True, encoding="utf-8", errors="replace", timeout=3600,
        )
    except subprocess.TimeoutExpired:
        log_lines.append(f"[error] {step} timed out (1 h limit)")
        _emit()
        st.warning(f"`{step}` timed out.")
        return False
    log_lines.extend((res.stdout or "").strip().splitlines()[-25:])
    log_lines.append(f"[exit {res.returncode}]")
    _emit()
    if res.returncode != 0:
        log_lines.extend((res.stderr or "").strip().splitlines()[-10:])
        _emit()
        st.warning(f"`{step}` exited with code {res.returncode}")
        return False
    return True


def _post_pipeline_steps(root: Path, run_dir: Path, log_lines: list[str],
                         log_box) -> None:
    """Generate the downstream report + publication figures after a run.

    Runs, in order: `visualize_results.py` (12 main figures),
    `prioritize_neural_tfs.py` (Track A/B shortlist + summary report) and
    `make_supp_go_figures.py` (4 GO supplementary figures + matrix CSV).
    Each step is gated on its inputs and runs as a subprocess so its CLI
    behaviour matches the documented commands exactly; output is teed into
    the same Streamlit log box as the pipeline itself.
    """
    st = _st()
    fig_out = root / "projects" / "NeuralTF" / "figures"
    top_csv = root / "projects" / "NeuralTF" / "results" / "top10_neural_tfs_prioritized.csv"

    st.markdown("**Generating the report + publication figures…**")

    _run_step(
        root, st, log_lines, log_box,
        ["projects/NeuralTF/scripts/visualize_results.py",
         "--run", str(run_dir), "--out", str(fig_out)],
        "visualize_results.py (12 main figures)",
        (run_dir / "rank.csv").exists(),
        "rank.csv not found in the run directory",
    )

    _run_step(
        root, st, log_lines, log_box,
        ["scripts/prioritize_neural_tfs.py", "--repo", str(root),
         "--rank", str(run_dir / "rank_neural.csv")],
        "prioritize_neural_tfs.py (Track A/B shortlist + report)",
        _downstream_readiness(root)["prioritization"]
        and (run_dir / "rank_neural.csv").exists(),
        "King 2024 mmc4/mmc5 xlsx, PlanMine parquet or rank_neural.csv missing",
    )

    _run_step(
        root, st, log_lines, log_box,
        ["projects/NeuralTF/scripts/make_supp_go_figures.py",
         "--run", str(run_dir),
         "--top-csv", str(top_csv),
         "--out", str(fig_out / "supplementary")],
        "make_supp_go_figures.py (4 GO supplementary figures)",
        top_csv.exists() and _downstream_readiness(root)["go_figures"],
        "top10_neural_tfs_prioritized.csv (from prioritization) or go.obo missing",
    )


def _data_rebuild_readiness(root: Path) -> dict:
    """Which raw -> processed build steps have their raw inputs available."""
    raw = root / "datasets" / "raw"
    king = raw / "Supplementary_Data_ King_2024"
    fincher_dge = (raw / "GSE111764_GEO_Fincher_atlas"
                   / "GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz")
    plass_tar = None
    if raw.exists():
        for p in raw.rglob("*.tar"):
            base = p.name.upper()
            if "RAW" in base or "GSE103633" in base:
                plass_tar = p
                break
    rosetta = raw / "smed_20140614.mapping.rosettastone.2020.txt"
    mmc4 = any(king.glob("*mmc4*.xlsx")) if king.exists() else False
    mmc7 = any(king.glob("*mmc7*.xlsx")) if king.exists() else False
    return {
        "fincher": fincher_dge.exists(),
        "plass": plass_tar is not None,
        "bridge": rosetta.exists() and mmc4,
        "king_atlas": mmc4 and mmc7,
        "planmine": (root / "projects" / "NeuralTF" / "runs" / "pipeline_run"
                     / "rank_neural.csv").exists(),
    }


def _rebuild_data(root: Path, log_lines: list[str], log_box) -> None:
    """Regenerate every processed data file from the raw downloads (gated)."""
    st = _st()
    reads = _data_rebuild_readiness(root)
    st.markdown("**Regenerating processed data files…**")

    _run_step(
        root, st, log_lines, log_box,
        ["scripts/convert_fincher.py"],
        "convert_fincher.py (fincher_subsample.h5ad)",
        reads["fincher"],
        "GSE111764 DGE .txt.gz not found under datasets/raw/GSE111764_GEO_Fincher_atlas/",
    )
    _run_step(
        root, st, log_lines, log_box,
        ["scripts/consolidate_plass.py"],
        "consolidate_plass.py (plass_v6.h5ad)",
        reads["plass"],
        "GSE103633_RAW.tar not found under datasets/raw/",
    )
    _run_step(
        root, st, log_lines, log_box,
        ["scripts/build_bridge.py"],
        "build_bridge.py (v4<->v6 bridge.csv)",
        reads["bridge"],
        "Rosetta Stone txt or King mmc4 xlsx missing",
    )
    _run_step(
        root, st, log_lines, log_box,
        ["scripts/build_king_atlas.py"],
        "build_king_atlas.py (king_atlas.tsv)",
        reads["king_atlas"],
        "King mmc4/mmc7 xlsx missing",
    )
    _run_step(
        root, st, log_lines, log_box,
        ["scripts/query_planmine.py", "--repo", str(root)],
        "query_planmine.py (PlanMine parquet + fasta, network)",
        reads["planmine"],
        "rank_neural.csv not found (run the pipeline first); "
        "this step queries the PlanMine web API",
    )


# ---------------------------------------------------------------------------
# Pages
# ---------------------------------------------------------------------------


def render_run_page() -> None:
    st = _st()
    st.subheader("Run the NeuralTF pipeline")

    status = _pipeline_status()
    ready = (
        status["fincher_h5ad"]
        and status["plass_h5ad"]
        and status["king_xlsx_ready"]
    )
    if not ready:
        st.warning("Required inputs are missing.")
        st.markdown(
            f"- `fincher_subsample.h5ad`: {'OK' if status['fincher_h5ad'] else 'MISSING'} "
            f"({status['fincher_path']})\n"
            f"- `plass_v6.h5ad`: {'OK' if status['plass_h5ad'] else 'MISSING'} "
            f"({status['plass_path']})\n"
            f"- King 2024 mmc4/mmc5/mmc6 xlsx: "
            f"{'OK' if status['king_xlsx_ready'] else 'MISSING'} "
            f"(`datasets/raw/Supplementary_Data_ King_2024/`)\n\n"
            "Prebuilt files come with the repo. Only the King 2024 xlsx tables "
            "are not committed - download mmc4-mmc7 from the Cell Reports paper "
            "supplementary. To rebuild the h5ads from raw GEO downloads:\n"
            "```\n"
            "python scripts/convert_fincher.py\n"
            "python scripts/consolidate_plass.py\n"
            "```"
        )
        st.caption(
            "Sources: GEO [GSE111764](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111764) "
            "(Fincher 2018) and [GSE103633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103633) "
            "(Plass 2018)."
        )
        return

    st.success("All inputs present. Ready to run.")
    ds = _downstream_readiness(_repo_root())
    st.caption(
        "Post-run generation (each step runs only if its inputs exist):\n"
        "- 12 main figures: OK\n"
        "- Track A/B shortlist + report: "
        f"{'OK' if ds['prioritization'] else 'MISSING (needs King 2024 mmc4/mmc5 xlsx + PlanMine parquet)'}\n"
        "- 4 GO supplementary figures: "
        f"{'OK' if ds['go_figures'] else 'MISSING (needs datasets/raw/go.obo)'}"
    )

    with st.form("pipeline_params"):
        col1, col2 = st.columns(2)
        with col1:
            subsample = st.number_input(
                "Cells per atlas (subsample)",
                min_value=0,
                max_value=40000,
                value=0,
                step=1000,
                help="0 = use the complete atlases (default); a value like "
                     "10000 speeds up development runs.",
            )
        with col2:
            out_dir = st.text_input(
                "Output directory (relative to repo root)",
                value="projects/NeuralTF/runs/pipeline_run",
                help="Where rank.csv, rank_neural.csv and evidence_cards.md are written.",
            )
        run_downstream = st.checkbox(
            "Also generate the report + all publication figures afterwards",
            value=True,
            help=("Runs prioritize_neural_tfs.py (Track A/B shortlist + candidate "
                  "summary report), visualize_results.py (9 figures) and "
                  "make_supp_go_figures.py (4 GO supplementary figures + matrix "
                  "CSV). Each step runs only when its inputs are available."),
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
        if run_downstream:
            _post_pipeline_steps(root, out_path, log_lines, log_box)
            for d in (
                root / "projects" / "NeuralTF" / "figures",
                root / "projects" / "NeuralTF" / "figures" / "supplementary",
                root / "projects" / "NeuralTF" / "results",
            ):
                if d.exists():
                    st.markdown(f"- `{d.relative_to(root)}`")
        st.button("View results →", on_click=lambda: st.session_state.update(page="Results"))
    except Exception as exc:  # noqa: BLE001
        progress.progress(0.0, text="Failed")
        st.error(f"Pipeline failed: {exc}")
        with st.expander("Log"):
            st.code("\n".join(log_lines[-200:]))

    _render_data_rebuild_expander()


def _render_data_rebuild_expander() -> None:
    """Regenerate every processed data file from the raw downloads in
    ``datasets/raw/`` (fincher, plass, bridge, king atlas, PlanMine)."""
    st = _st()
    root = _repo_root()
    with st.expander("Regenerate processed data files (from raw downloads)",
                     expanded=False):
        reads = _data_rebuild_readiness(root)
        st.caption(
            "Rebuilds every processed file from the raw sources in "
            "`datasets/raw/`: fincher/plass h5ads, bridge.csv, king_atlas.tsv "
            "and the PlanMine parquet+fasta (network). Each step runs only if "
            "its raw input is present. All files are generated locally — "
            "use this after downloading the raw sources to their documented "
            "locations (or run `python scripts/generate_all.py`)."
        )
        st.markdown(
            "\n".join(
                f"- {k}: {'OK' if v else 'MISSING'}"
                for k, v in reads.items()
            )
        )
        if st.button("Rebuild processed data files now", type="secondary"):
            log_box = st.empty()
            _rebuild_data(root, [], log_box)


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

    # Tabbed view of the two rankings + cards + run metadata + visualizations
    tab_rank, tab_neural, tab_cards, tab_viz, tab_supp, tab_meta, tab_dirichlet, tab_fstf = st.tabs(
        ["All candidates", "Neural-filtered", "Evidence cards",
         "Visualizations", "GO supplementary", "Run metadata",
         "Dirichlet sensitivity", "FSTF rankings"]
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

    with tab_viz:
        if rank_csv.exists():
            _render_visualizations(pd.read_csv(rank_csv), run_dir=run)
        else:
            st.warning("rank.csv not found in this run.")

    with tab_supp:
        _render_go_supplementary(_repo_root())

    with tab_meta:
        if json_path.exists():
            try:
                meta = json.loads(json_path.read_text(encoding="utf-8"))
                st.json(meta)
            except Exception as exc:  # noqa: BLE001
                st.warning(f"Could not parse {json_path.name}: {exc}")
        else:
            st.info("No pipeline_results.json in this run.")

    with tab_dirichlet:
        _render_dirichlet_sensitivity(_repo_root())

    with tab_fstf:
        _render_fstf_rankings(_repo_root())


# ---------------------------------------------------------------------------
# Visualizations tab — reuses the figure builders from
# projects/NeuralTF/scripts/visualize_results.py without re-saving PNGs
# to disk. Streamlit can render a returned matplotlib Figure directly via
# st.pyplot(fig). The script-side wrappers (the `fig_*` functions) save
# the same figure to disk; here we use the `make_*` builders instead.
# ---------------------------------------------------------------------------


def _import_visualize_fixed():
    """Import visualize_fixed.py from projects/NeuralTF/scripts/."""
    import importlib.util
    p = _repo_root() / "projects" / "NeuralTF" / "scripts" / "visualize_fixed.py"
    if not p.exists():
        return None
    spec = importlib.util.spec_from_file_location("visualize_fixed", p)
    if spec is None or spec.loader is None:
        return None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _call_make(func, *args, **kwargs):
    """Call a make_* function and return the figure without saving."""
    import tempfile
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with tempfile.TemporaryDirectory() as tmpdir:
        out_path = Path(tmpdir)
        # Call the function with out_path
        func(*args, out_path=out_path, **kwargs)
        # Get the figure that was created
        figs = [plt.figure(n) for n in plt.get_fignums()]
        if figs:
            fig = figs[-1]  # Get the last created figure
            return fig
    return None


def _render_visualizations(df, *, run_dir) -> None:
    st = _st()
    mod = _import_visualize_fixed()
    if mod is None:
        st.warning(
            "Could not import `projects/NeuralTF/scripts/visualize_fixed.py`."
            " Make sure it exists, or run: `python projects/NeuralTF/scripts/visualize_fixed.py` as a script."
        )
        return

    import matplotlib
    import pandas as pd
    previous_backend = matplotlib.get_backend()
    try:
        matplotlib.use("Agg", force=True)
        c1, c2, c3 = st.columns(3)
        c1.metric("Candidates", len(df))
        if "tier" in df.columns:
            c2.metric("HIGH tier",
                      int((df["tier"].str.lower() == "high").sum()))
            c3.metric("RNAi-validated",
                      int((df.get("rnai", pd.Series([])) > 0).sum())
                      if "rnai" in df.columns else 0)
        st.caption(
            "*These figures are derived from this run's `rank.csv` (and the "
            "GO-term matrix in `figures/supplementary/` for the GO panel) — "
            "no hardcoded meta-data, no external lookup tables. Save PNGs via "
            "*`python projects/NeuralTF/scripts/visualize_fixed.py`*."
        )

        neural_df = None
        neural_csv = run_dir / "rank_neural.csv"
        if neural_csv.exists():
            neural_df = pd.read_csv(neural_csv)

        builders = [
            ("Candidate summary (tiers / proof / score / coverage)",
             mod.make_candidate_summary, [df]),
            ("Top-10 dual track (Track A vs Track B)",
             mod.make_top10_dual_track, [neural_df]),
            ("Evidence matrix (top 30)",
             mod.make_evidence_heatmap, [df], {"n": 30}),
            ("Candidate funnel (scored → neural → final)",
             mod.make_candidate_funnel, [df, neural_df]),
            ("Evidence composition (top 15)",
             mod.make_evidence_composition, [df], {"n": 15}),
            ("Stream ablation (rank sensitivity)",
             mod.make_stream_ablation, [df]),
            ("Top-10 radar fingerprints",
             mod.make_top10_radar, [neural_df]),
            ("GO dot plot (top-10 terms)",
             mod.make_go_dotplot, [neural_df]),
            ("Score distributions (all streams)",
             mod.make_score_distributions, [df]),
            ("Integrated vs composite (bonuses)",
             mod.make_integrated_vs_composite, [neural_df]),
            ("Proof-status violin (score distributions)",
             mod.make_proof_status_violin,
             [neural_df if neural_df is not None else df]),
            ("Weight sensitivity (Top-10 rank bands)",
             mod.make_weight_sensitivity, [neural_df]),
            ("Integrated score vs neural filter (ECDF)",
             mod.make_integrated_vs_neural_filter, [df, neural_df]),
        ]
        left, right = st.columns(2)
        try:
            for i, entry in enumerate(builders):
                col = left if i % 2 == 0 else right
                title = entry[0]
                fn = entry[1]
                args = entry[2]
                kwargs = entry[3] if len(entry) > 3 else {}
                with col:
                    st.markdown(f"**{title}**")
                    try:
                        fig = _call_make(fn, *args, **kwargs)
                    except Exception as exc:
                        st.warning(f"Could not build `{title}`: {exc}")
                        continue
                    if fig is None:
                        st.caption("_(skipped — required columns missing)_")
                    else:
                        st.pyplot(fig)
                        plt_close(fig)
        finally:
            matplotlib.use(previous_backend, force=True)
    finally:
        matplotlib.use(previous_backend, force=True)


def _render_dirichlet_sensitivity(root: Path) -> None:
    """Render the Dirichlet sensitivity analysis outputs.

    Shows centered/uniform/all249 CSVs and method comparison figures.
    """
    st = _st()
    import pandas as pd
    res = root / "projects" / "NeuralTF" / "results"
    fig_dir = root / "projects" / "NeuralTF" / "figures"

    st.subheader("Dirichlet sensitivity analysis")
    st.caption(
        "Weight-sensitivity analysis using Dirichlet sampling. "
        "Three methods: fixed-weight baseline, centered Dirichlet (k=40), "
        "and uniform Dirichlet (alpha=1)."
    )

    # --- Centered Dirichlet ------------------------------------------------
    with st.expander("Centered Dirichlet (k=40, informative prior)", expanded=False):
        centered_csv = res / "dirichlet_top10_prioritized.csv"
        overall_csv = res / "dirichlet_overall_top10.csv"
        report_md = res / "dirichlet_candidate_summary_report.md"

        if centered_csv.exists():
            df = pd.read_csv(centered_csv)
            st.markdown("**Track-based top-10** (5A + 5B)")
            st.dataframe(df, use_container_width=True, hide_index=True)
        else:
            st.info("dirichlet_top10_prioritized.csv not found. Run `python projects/NeuralTF/scripts/dirichlet_prioritize.py`.")

        if overall_csv.exists():
            df_overall = pd.read_csv(overall_csv)
            st.markdown("**Overall top-10 by Dirichlet median score**")
            st.dataframe(df_overall, use_container_width=True, hide_index=True)

        if report_md.exists():
            with st.expander("Candidate summary report"):
                st.markdown(report_md.read_text(encoding="utf-8"))

        # Centered figures
        centered_figs = [
            "fig_dirichlet_trackA_top5.png",
            "fig_dirichlet_trackB_top5.png",
            "fig_dirichlet_scatter.png",
            "fig_dirichlet_combined.png",
            "fig_dirichlet_score_shift.png",
        ]
        for fname in centered_figs:
            fpath = fig_dir / fname
            if fpath.exists():
                st.image(str(fpath), use_container_width=True, caption=fname.replace(".png", "").replace("_", " "))

    # --- Uniform Dirichlet -------------------------------------------------
    with st.expander("Uniform Dirichlet (alpha=1, non-informative prior)", expanded=False):
        unif_csv = res / "dirichlet_uniform_top10.csv"
        unif_full = res / "dirichlet_uniform_full_rank.csv"
        unif_summary = res / "dirichlet_uniform_summary.txt"

        if unif_csv.exists():
            df_unif = pd.read_csv(unif_csv)
            st.markdown("**Track-based top-10** (5A + 5B, uniform prior)")
            st.dataframe(df_unif, use_container_width=True, hide_index=True)
        else:
            st.info("dirichlet_uniform_top10.csv not found. Run `python projects/NeuralTF/scripts/dirichlet_uniform.py`.")

        if unif_full.exists():
            df_full = pd.read_csv(unif_full)
            st.caption(f"Full rank: {len(df_full)} candidates with uniform scores")

        if unif_summary.exists():
            with st.expander("Summary"):
                st.code(unif_summary.read_text(encoding="utf-8"))

        # Uniform figures
        unif_figs = [
            "fig_dirichlet_uniform_trackA_top5.png",
            "fig_dirichlet_uniform_trackB_top5.png",
            "fig_dirichlet_uniform_scatter.png",
            "fig_dirichlet_uniform_combined.png",
            "fig_dirichlet_uniform_score_shift.png",
        ]
        for fname in unif_figs:
            fpath = fig_dir / fname
            if fpath.exists():
                st.image(str(fpath), use_container_width=True, caption=fname.replace(".png", "").replace("_", " "))

    # --- 99 vs 249 comparison ----------------------------------------------
    with st.expander("99-neural vs 249-wide comparison", expanded=False):
        all249_csv = res / "dirichlet_uniform_all249_top10.csv"
        all249_summary = res / "dirichlet_uniform_all249_summary.txt"
        fig_99vs249 = fig_dir / "fig_dirichlet_99vs249.png"

        if all249_csv.exists():
            df_249 = pd.read_csv(all249_csv)
            st.markdown("**249-wide top-10** (uniform Dirichlet, no neural filter)")
            st.dataframe(df_249, use_container_width=True, hide_index=True)
        else:
            st.info("dirichlet_uniform_all249_top10.csv not found. Run `python projects/NeuralTF/scripts/dirichlet_uniform_all249.py`.")

        if all249_summary.exists():
            st.code(all249_summary.read_text(encoding="utf-8"))

        if fig_99vs249.exists():
            st.image(str(fig_99vs249), use_container_width=True, caption="99-neural vs 249-wide rank shift + score comparison")

    # --- 3-way method comparison -------------------------------------------
    with st.expander("3-way method comparison figures", expanded=False):
        comp_figs = [
            "fig_dirichlet_3way_comparison.png",
            "fig_dirichlet_uniform_vs_centered.png",
            "fig_dirichlet_score_density.png",
            "fig_dirichlet_rank_correlation.png",
            "fig_dirichlet_score_volatility.png",
            "fig_dirichlet_method_summary.png",
        ]
        found_any = False
        for fname in comp_figs:
            fpath = fig_dir / fname
            if fpath.exists():
                st.image(str(fpath), use_container_width=True, caption=fname.replace(".png", "").replace("_", " "))
                found_any = True
        if not found_any:
            st.info("Method comparison figures not found. Run `python projects/NeuralTF/scripts/dirichlet_method_comparison.py`.")


def _render_fstf_rankings(root: Path) -> None:
    """Render the three FSTF ranking CSVs (19/43/74 scope)."""
    st = _st()
    import pandas as pd
    res = root / "projects" / "NeuralTF" / "results"

    st.subheader("FSTF (Planarian Stem Cell TF) rankings")
    st.caption(
        "Three scope levels from King 2024 mmc4 TF catalog (FSTF? = yes). "
        "All sorted by composite score (descending)."
    )

    scopes = [
        ("fstf_ranked_19_neural.csv", "19 FSTFs — neural-filtered",
         "FSTFs with King neural signal or RNAi evidence. Most relevant to the thesis question."),
        ("fstf_ranked_43_all.csv", "43 FSTFs — all candidates",
         "FSTFs that passed the expression filter (p <= 0.05)."),
        ("fstf_ranked_74_catalog.csv", "74 FSTFs — full catalog",
         "All FSTFs from King mmc4 TF sheet, regardless of expression filter."),
    ]

    for fname, title, desc in scopes:
        fpath = res / fname
        with st.expander(title, expanded=(fname == "fstf_ranked_19_neural.csv")):
            st.caption(desc)
            if fpath.exists():
                df = pd.read_csv(fpath)
                st.dataframe(df, use_container_width=True, hide_index=True)
            else:
                st.info(f"`{fname}` not found. Run `python projects/NeuralTF/scripts/export_fstf_ranked.py`.")


def _render_go_supplementary(root: Path) -> None:
    """Show the 4 GO supplementary figures + matrix from figures/supplementary.

    The figures are PNG artifacts generated by `make_supp_go_figures.py`
    (run automatically after a pipeline run, or via
    `python scripts/generate_all.py`).
    """
    st = _st()
    supp = root / "projects" / "NeuralTF" / "figures" / "supplementary"
    pngs = sorted(supp.glob("fig_s*.png")) if supp.exists() else []
    if not pngs:
        st.info(
            "GO supplementary figures are not generated yet. Run the pipeline "
            "(**Run** page, the 'Also generate the report + publication "
            "figures' option) or: `python scripts/generate_all.py`."
        )
        return
    st.caption(
        "GO gene-term map (S1), Top-10 dot matrix (S2), top/shared-term "
        "statistics (S3) and neural-only focus (S4), produced from the "
        "PlanMine annotations + go.obo."
    )
    for png in pngs:
        st.image(str(png), use_container_width=True, caption=png.stem)
    mat = supp / "go_gene_term_matrix_reduced.csv"
    if mat.exists():
        import pandas as pd
        df = pd.read_csv(mat)
        st.caption(f"Reduced gene × term matrix: {df.shape[0]} genes × "
                   f"{df.shape[1] - 1} terms")
        with st.expander("Show matrix"):
            st.dataframe(df, use_container_width=True)


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


def render_prioritization_page() -> None:
    """Show the dual-track shortlist from all three methods as tabs."""
    st = _st()
    import pandas as pd
    root = _repo_root()
    res = root / "projects" / "NeuralTF" / "results"

    st.subheader("Prioritized neural-fate TF shortlist")

    tab_fixed, tab_centered, tab_uniform = st.tabs([
        "Fixed-weight (baseline)",
        "Dirichlet-centered (k=40)",
        "Dirichlet-uniform (alpha=1)",
    ])

    # --- Fixed-weight tab ---------------------------------------------------
    with tab_fixed:
        csv_fixed = res / "top10_neural_tfs_prioritized.csv"
        if csv_fixed.exists():
            df = pd.read_csv(csv_fixed)
            st.caption(f"{len(df)} TFs (5A + 5B) — fixed-weight baseline")
            ta = df[df["track"] == "A"]
            tb = df[df["track"] == "B"]
            st.markdown("**Track A — RNAi-validated**")
            st.dataframe(ta, use_container_width=True, hide_index=True)
            st.markdown("**Track B — novel candidates**")
            st.dataframe(tb, use_container_width=True, hide_index=True)
        else:
            st.info("Not generated yet. Run `python scripts/prioritize_neural_tfs.py`.")

    # --- Dirichlet-centered tab ---------------------------------------------
    with tab_centered:
        csv_centered = res / "dirichlet_top10_prioritized.csv"
        overall_centered = res / "dirichlet_overall_top10.csv"
        report_md = res / "dirichlet_candidate_summary_report.md"

        if csv_centered.exists():
            df = pd.read_csv(csv_centered)
            st.caption(f"{len(df)} TFs (5A + 5B) — Dirichlet-centered (k=40)")
            ta = df[df["track"] == "A"]
            tb = df[df["track"] == "B"]
            st.markdown("**Track A — RNAi-validated**")
            st.dataframe(ta, use_container_width=True, hide_index=True)
            st.markdown("**Track B — novel candidates**")
            st.dataframe(tb, use_container_width=True, hide_index=True)
        else:
            st.info("Not generated yet. Run `python projects/NeuralTF/scripts/dirichlet_prioritize.py`.")

        if overall_centered.exists():
            with st.expander("Overall top-10 by Dirichlet median score"):
                st.dataframe(pd.read_csv(overall_centered), use_container_width=True, hide_index=True)

        if report_md.exists():
            with st.expander("Candidate summary report"):
                st.markdown(report_md.read_text(encoding="utf-8"))

    # --- Dirichlet-uniform tab ----------------------------------------------
    with tab_uniform:
        csv_uniform = res / "dirichlet_uniform_top10.csv"
        overall_uniform = res / "dirichlet_uniform_overall_top10.csv"
        summary_uniform = res / "dirichlet_uniform_summary.txt"

        if csv_uniform.exists():
            df = pd.read_csv(csv_uniform)
            st.caption(f"{len(df)} TFs (5A + 5B) — Dirichlet-uniform (alpha=1)")
            ta = df[df["track"] == "A"]
            tb = df[df["track"] == "B"]
            st.markdown("**Track A — RNAi-validated**")
            st.dataframe(ta, use_container_width=True, hide_index=True)
            st.markdown("**Track B — novel candidates**")
            st.dataframe(tb, use_container_width=True, hide_index=True)
        else:
            st.info("Not generated yet. Run `python projects/NeuralTF/scripts/dirichlet_uniform.py`.")

        if overall_uniform.exists():
            with st.expander("Overall top-10 by uniform median score"):
                st.dataframe(pd.read_csv(overall_uniform), use_container_width=True, hide_index=True)

        if summary_uniform.exists():
            with st.expander("Summary"):
                st.code(summary_uniform.read_text(encoding="utf-8"))


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
        "discovery (3 atlases, 7 evidence streams)."
    )
    page = st.sidebar.selectbox("Page", ["Run", "Results", "Prioritization", "Assistant"])
    if page == "Run":
        render_run_page()
    elif page == "Results":
        render_results_page()
    elif page == "Prioritization":
        render_prioritization_page()
    else:
        render_assistant_page()


if __name__ == "__main__":  # pragma: no cover — exercised via `streamlit run`
    if "--help" in sys.argv:
        print("Streamlit app. Launch with: bioforge ui  or  streamlit run src/bioforge/ui/app.py")
        sys.exit(0)
    main()
