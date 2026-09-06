"""Optional FASTQ -> matrix orchestration (stub).

The actual FASTQ->matrix recipe requires installing `kb-python` (kallisto |
  bustools) for 10x-style or `salmon alevin` for non-UMI bulk via the
  optional `[fastq]` extra. This module declares the recipe API and a
  minimal hook so workflows can opt in via `recipe: fastq_to_10x`.

The first implementation is intentionally a stub: it raises a clear error
unless the recipe is implemented in this env. The stub keeps the import
lightweight so the default BioForge container doesn't pull in heavy
aligners.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class FastqRecipe:
    name: str
    description: str

    def is_implemented(self) -> bool:
        return False  # explicit: drug of the day not wired up yet

    def execute(self, fastq_dir: Path, output_dir: Path,
                transcriptome_fasta: Optional[Path] = None) -> Path:
        raise NotImplementedError(
            f"FASTQ recipe '{self.name}' is not implemented in this BioForge "
            "install. Install the optional [fastq] extra (kb-python OR "
            "salmon) and see docs/recipes/ for guidance."
        )


FASTQ_TO_10X = FastqRecipe(
    name="fastq_to_10x",
    description="Convert paired-end FASTQ into a 10x-style mtx directory via kb-python."
)


def available_recipes() -> list[str]:
    return [FASTQ_TO_10X.name]
