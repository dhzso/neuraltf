"""Master script to regenerate all publication figures.

Usage:
    python projects/NeuralTF/scripts/generate_publication_figures.py
    python projects/NeuralTF/scripts/generate_publication_figures.py --figure 1 3

Outputs:
    Main figures:    projects/NeuralTF/figures/Fig{1-6}_*.png + .pdf
    Supplementary:   projects/NeuralTF/figures/supplementary/fig_s{1-4}_*.png
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from pathlib import Path
import traceback

REPO = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = Path(__file__).resolve().parent
FIGURES_DIR = SCRIPTS_DIR / "figures"
SUPP_DIR = REPO / "projects" / "NeuralTF" / "figures" / "supplementary"
SUPP_DIR.mkdir(parents=True, exist_ok=True)

FIGURE_MODULES = {
    1: ("Fig 1 — Candidate landscape",          FIGURES_DIR / "Fig1_candidate_landscape.py"),
    2: ("Fig 2 — Evidence architecture",         FIGURES_DIR / "Fig2_evidence_architecture.py"),
    3: ("Fig 3 — Ranking robustness",            FIGURES_DIR / "Fig3_ranking_robustness.py"),
    4: ("Fig 4 — Stream sensitivity",            FIGURES_DIR / "Fig4_stream_sensitivity.py"),
    5: ("Fig 5 — Neural filtering scope",        FIGURES_DIR / "Fig5_neural_vs_full.py"),
    6: ("Fig 6 — Prioritized candidate atlas",   FIGURES_DIR / "Fig6_prioritized_candidate_atlas.py"),
}


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    parser = argparse.ArgumentParser(description="Regenerate publication figures")
    parser.add_argument("--figure", nargs="*", type=int,
                        help="Specific figures to generate (default: all)")
    args = parser.parse_args()

    t0 = time.time()
    generated = []
    errors = []

    fig_nums = args.figure if args.figure else sorted(FIGURE_MODULES.keys())

    for num in fig_nums:
        if num not in FIGURE_MODULES:
            print(f"  [SKIP] Figure {num} not available")
            continue
        label, path = FIGURE_MODULES[num]
        print(f"\n{'='*60}")
        print(f"  Generating {label}")
        print(f"{'='*60}")
        try:
            mod = _load_module(path, f"fig{num}")
            mod.build()
            generated.append(num)
            print(f"  OK Figure {num} complete")
        except Exception:
            print(f"  FAILED Figure {num}:")
            traceback.print_exc()
            errors.append(num)

    # Supplementary GO figures
    print(f"\n{'='*60}")
    print(f"  Generating supplementary GO figures")
    print(f"{'='*60}")
    supp_script = SCRIPTS_DIR / "make_supp_go_figures.py"
    if supp_script.exists():
        try:
            mod = _load_module(supp_script, "make_supp_go_figures")
            if hasattr(mod, "build_all"):
                mod.build_all()
            elif hasattr(mod, "main"):
                mod.main()
            print(f"  OK Supplementary figures complete")
        except Exception:
            print(f"  FAILED Supplementary figures:")
            traceback.print_exc()
    else:
        print(f"  SKIP make_supp_go_figures.py not found")

    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"  SUMMARY")
    print(f"{'='*60}")
    print(f"  Generated: {len(generated)} main figures")
    print(f"  Errors:    {len(errors)}")
    print(f"  Time:      {elapsed:.1f}s")

    if errors:
        print(f"\n  Failed figures: {errors}")
        sys.exit(1)


if __name__ == "__main__":
    main()
