#!/usr/bin/env python
"""Post-pipeline downstream orchestration.

Runs all scoring, Dirichlet, ANANSE, supplementary table, and figure
generation steps in dependency order after the main pipeline completes.

Usage:
    .venv\Scripts\python.exe scripts\run_downstream.py [--force]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PYTHON = ROOT / ".venv" / "Scripts" / "python.exe"
if not PYTHON.exists():
    PYTHON = Path(sys.executable)


def run(script_parts: list, label: str, force: bool = False, outputs: list[Path] | None = None) -> bool:
    """Run a script; skip if all outputs already exist (unless --force)."""
    if not force and outputs:
        if all(p.exists() and p.stat().st_size > 0 for p in outputs):
            print(f"  [SKIP] {label} — all outputs present")
            return True

    script = ROOT.joinpath(*script_parts)
    if not script.exists():
        print(f"  [WARN] {label} — script not found: {script}")
        return False

    print(f"\n>>> {label}")
    t0 = time.time()
    result = subprocess.run(
        [str(PYTHON), str(script)],
        cwd=str(ROOT),
        capture_output=False,
    )
    elapsed = time.time() - t0
    if result.returncode == 0:
        print(f"    [OK] {label} ({elapsed:.0f}s)")
        return True
    else:
        print(f"    [FAIL] {label} exited {result.returncode} ({elapsed:.0f}s)")
        return False


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--force", action="store_true", help="Re-run even if outputs exist")
    args = p.parse_args()
    force = args.force

    RUN = ROOT / "projects" / "NeuralTF" / "runs" / "pipeline_run"
    RES = ROOT / "projects" / "NeuralTF" / "results"
    FIG = ROOT / "projects" / "NeuralTF" / "figures"

    # Verify pipeline outputs exist
    rank_csv = RUN / "rank.csv"
    if not rank_csv.exists():
        print(f"[ERROR] rank.csv not found at {rank_csv}. Run scripts/run.py first.")
        return 1

    import pandas as pd
    rank = pd.read_csv(rank_csv)
    n_cand = len(rank)
    print(f"[INFO] rank.csv: {n_cand} candidates")
    if "perez_lineage" in rank.columns:
        n_perez = rank["perez_lineage"].notna().sum()
        print(f"[INFO] perez_lineage column present: {n_perez} non-null values")
    else:
        print("[WARN] perez_lineage column MISSING from rank.csv — pipeline may not have included Perez stream")

    errors = []

    steps = [
        # (label, script_parts, output_files)
        ("Dirichlet Centered (n=99)",
         ["projects", "NeuralTF", "scripts", "dirichlet_prioritize.py"],
         [RES / "dirichlet_centered_full_rank.csv"]),

        ("Dirichlet Uniform (n=99)",
         ["projects", "NeuralTF", "scripts", "dirichlet_uniform.py"],
         [RES / "dirichlet_uniform_full_rank.csv"]),

        ("Dirichlet Uniform All-249",
         ["projects", "NeuralTF", "scripts", "dirichlet_uniform_all249.py"],
         [RES / "dirichlet_uniform_all249_full_rank.csv"]),

        ("Dirichlet Centered All-249",
         ["projects", "NeuralTF", "scripts", "dirichlet_centered_all249.py"],
         [RES / "dirichlet_centered_all249_full_rank.csv"]),

        ("Export ranked FSTF",
         ["projects", "NeuralTF", "scripts", "export_fstf_ranked.py"],
         [RES / "fstf_ranked_19_neural.csv"]),

        ("ANANSE full scan",
         ["projects", "NeuralTF", "scripts", "ananse_full_scan.py"],
         [RES / "ananse_network_full.csv"]),

        ("Supplementary tables",
         ["projects", "NeuralTF", "scripts", "create_supplementary_tables.py"],
         [RES / "supplementary_table_S1_method_comparison.csv"]),

        ("Publication figures",
         ["projects", "NeuralTF", "scripts", "generate_publication_figures.py"],
         [FIG / "01_stream_coverage_249.png"]),
    ]

    for label, script_parts, outputs in steps:
        ok = run(script_parts, label, force=force, outputs=[Path(o) for o in outputs])
        if not ok:
            errors.append(label)

    print("\n" + "=" * 60)
    if errors:
        print(f"[DONE WITH ERRORS] {len(errors)} steps failed:")
        for e in errors:
            print(f"  - {e}")
        return 1
    else:
        print(f"[DONE] All {len(steps)} downstream steps completed successfully.")
        return 0


if __name__ == "__main__":
    sys.exit(main())
