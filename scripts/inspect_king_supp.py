"""Inspect King 2024 supplementary xlsx files to learn the real schema.

Usage (inside container):
    python /workspace/scripts/inspect_king_supp.py
"""
from __future__ import annotations

import os
import sys

import pandas as pd

BASE = "datasets/raw/Supplementary_Data_ King_2024"


def inspect_excel(path: str, max_sheets: int = 4, n_rows: int = 6) -> None:
    print("=" * 78)
    print("FILE:", os.path.basename(path))
    xl = pd.ExcelFile(path)
    print("sheets:", xl.sheet_names)
    for sh in xl.sheet_names[:max_sheets]:
        df = pd.read_excel(path, sheet_name=sh, header=None, nrows=n_rows)
        print(f"--- sheet: {sh!r}  shape={df.shape}")
        with pd.option_context("display.max_columns", 10, "display.width", 200):
            print(df.to_string())
        print()


def main() -> int:
    files = sorted(f for f in os.listdir(BASE) if f.endswith(".xlsx"))
    for f in files:
        try:
            inspect_excel(os.path.join(BASE, f))
        except Exception as exc:
            print("ERROR reading", f, ":", exc)
    return 0


if __name__ == "__main__":
    sys.exit(main())
