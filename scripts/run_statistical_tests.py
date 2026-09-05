#!/usr/bin/env python
"""Run all statistical tests for the NeuralTF pipeline.

Runs all 14 statistical tests in dependency order, enforcing
prerequisites (the Dirichlet draw matrices several tests consume) and a
UTF-8 child environment (Windows cp1252 consoles cannot encode the
arrow/Greek characters the stats scripts print).
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
STATS_DIR = REPO / "scripts" / "stats"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"

CHILD_ENV = {**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"}

TESTS = [
    ("permutation_test_full.py", []),
    ("bootstrap_confidence.py", []),
    ("overlap_significance.py", []),
    ("precision_recall.py", []),
    ("negative_controls.py", []),
    ("effect_sizes.py", []),
    ("leave_one_atlas_out.py", []),
    ("meta_analytic_pvalue.py", []),
    ("power_analysis.py", []),
    ("mann_whitney_top10.py", []),
    ("calibration.py", []),
    ("brier_score.py", []),
    ("cross_method_correction.py", []),
    ("score_shuffling_permutation.py", ["--n-perm", "200"]),
]


def check_prerequisites() -> list[str]:
    """Fail fast with actionable messages when inputs are missing."""
    missing = []
    rank_csv = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run" / "rank.csv"
    if not rank_csv.exists():
        missing.append(f"rank.csv not found at {rank_csv} - run scripts/run.py first")
    # bootstrap_confidence / power_analysis / overlap_significance consume
    # the Dirichlet draw/full-rank artifacts; running stats without them
    # silently degrades bootstrap to its parametric fallback.
    for name in (
        "dirichlet_centered_draw_scores.csv",
        "dirichlet_uniform_draw_scores.csv",
        "dirichlet_centered_full_rank.csv",
        "dirichlet_uniform_full_rank.csv",
    ):
        p = RESULTS_DIR / name
        if not p.exists():
            missing.append(
                f"{p.name} missing - run projects/NeuralTF/scripts/dirichlet_centered.py "
                f"and dirichlet_uniform.py first"
            )
    return missing


def main() -> int:
    print("=== Running All Statistical Tests ===", flush=True)

    missing = check_prerequisites()
    if missing:
        for m in missing:
            print(f"  [ERROR] {m}", flush=True)
        return 1

    # The 3-atlas production-scale permutation (re-clusters 143k cells and
    # re-runs Wilcoxon DE x 30 permutations) takes ~2.5-3 h; the driver
    # default timeout is 1 h, so the permutation gets its own budget.
    TIMEOUTS = {"permutation_test_full.py": 4 * 3600}

    results = []
    for script_name, extra_args in TESTS:
        script_path = STATS_DIR / script_name
        if not script_path.exists():
            results.append((script_name, "skipped (script not found)"))
            continue
        print(f"\n>>> {script_name}", flush=True)
        try:
            res = subprocess.run(
                [sys.executable, str(script_path), *extra_args],
                cwd=REPO,
                timeout=TIMEOUTS.get(script_name, 3600),
                env=CHILD_ENV,
            )
            results.append((script_name, "OK" if res.returncode == 0 else f"FAILED (exit {res.returncode})"))
        except subprocess.TimeoutExpired:
            results.append((script_name, "FAILED (timeout)"))

    print("\n=== Summary ===", flush=True)
    for name, status in results:
        print(f"  {name:<40} {status}", flush=True)

    failed = [r for r in results if r[1].startswith("FAILED")]
    if failed:
        print(f"\n{len(failed)} test(s) failed.", flush=True)
        return 1
    print("\nAll tests completed successfully.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
