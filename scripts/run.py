#!/usr/bin/env python
"""Run the NeuralTF pipeline end-to-end.

Thin wrapper around the CLI; useful when you want to invoke the pipeline
from a script (e.g. inside an editor or a notebook) without going through
the `bioforge` console script.

Usage:
    python scripts/run.py [--subsample 0] [--out DIR]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from bioforge.projects.neuraltf.pipeline import NeuralTFPipeline  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--subsample",
        type=int,
        default=0,
        help="Cells per atlas (default 0 = keep the complete atlases; "
             "e.g. 10000 for a fast development run)",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output directory (default projects/NeuralTF/runs/pipeline_run)",
    )
    args = p.parse_args()

    pipe = NeuralTFPipeline(
        out_dir=args.out,
        subsample=args.subsample if args.subsample > 0 else None,
    )
    pipe.run()


if __name__ == "__main__":
    main()
