#!/usr/bin/env python
"""Run all statistical tests for the NeuralTF pipeline.

This script runs all 14 statistical tests in the correct order,
ensuring dependencies are met.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
STATS_DIR = REPO / "scripts" / "stats"

TESTS = [
    ("permutation_test_full.py", ["--n-perm", "100", "--subsample", "2000"]),
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


def main() -> int:
    print("=== Running All Statistical Tests ===", flush=True)
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
                timeout=3600,
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
