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
    ("permutation_test_full.py", ["--n-perm", "1000", "--subsample", "2000"]),
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
    ("score_shuffling_permutation.py", []),
]


def main() -> int:
    print("=== Running All Statistical Tests ===")
    results = []
    for script_name, extra_args in TESTS:
        script_path = STATS_DIR / script_name
        if not script_path.exists():
            results.append((script_name, "skipped (script not found)"))
            continue
        print(f"\n>>> {script_name}")
        try:
            res = subprocess.run(
                [sys.executable, str(script_path), *extra_args],
                cwd=REPO,
                timeout=3600,
            )
            results.append((script_name, "OK" if res.returncode == 0 else f"FAILED (exit {res.returncode})"))
        except subprocess.TimeoutExpired:
            results.append((script_name, "FAILED (timeout)"))
    
    print("\n=== Summary ===")
    for name, status in results:
        print(f"  {name:<40} {status}")
    
    failed = [r for r in results if r[1].startswith("FAILED")]
    if failed:
        print(f"\n{len(failed)} test(s) failed.")
        return 1
    print("\nAll tests completed successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
