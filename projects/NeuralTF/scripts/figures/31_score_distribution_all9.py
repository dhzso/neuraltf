"""Legacy wrapper for 31_score_distribution_all11.py (backwards compatibility)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib.util

spec = importlib.util.spec_from_file_location("fig31_11", Path(__file__).parent / "31_score_distribution_all11.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
build = mod.build

if __name__ == "__main__":
    build()
