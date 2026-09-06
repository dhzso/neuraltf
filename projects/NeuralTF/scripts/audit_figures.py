#!/usr/bin/env python
"""Audit all publication figures for non-emptiness and completeness.

A figure is EMPTY when its PNG is dominated by background pixels — i.e. an
axes-only plot, a placeholder text panel, or a file written with no data.
This tool:

  1. walks projects/NeuralTF/figures (main + supplementary)
  2. flags files < MIN_BYTES (truncated writes) and files whose
     non-background ink coverage is below threshold (empty axes)
  3. verifies the expected 33 numbered figures + 7 supplementary figures
     all exist
  4. writes figures/figure_completeness_report.json

Usage:
    python projects/NeuralTF/scripts/audit_figures.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
FIG = REPO / "projects" / "NeuralTF" / "figures"
SUP = FIG / "supplementary"

MIN_BYTES = 15_000          # below this a 300-dpi PNG is almost surely axes-only
INK_COVERAGE_MIN = 0.005    # <0.5% non-background pixels => effectively empty

EXPECTED_MAIN = [f"{i:02d}" for i in range(1, 34)]
EXPECTED_SUPP = [
    "fig_s1_go_gene_term_map", "fig_s2_go_top10_dotmatrix",
    "fig_s3_go_top_terms", "fig_s4_go_neural_focus",
    "fig_s5_go_heatmap_neural", "fig_s6_top10_go_profiles",
    "fig_s7_go_namespace_track",
]


def _ink_coverage(png_path: Path) -> float | None:
    """Fraction of pixels that differ from the (dominant) background color."""
    try:
        from PIL import Image
    except ImportError:
        return None
    try:
        with Image.open(png_path) as im:
            im = im.convert("RGB")
            im.thumbnail((400, 400))
            arr = np.asarray(im, dtype=np.int16)
    except Exception:
        return None
    # background = most common color
    flat = arr.reshape(-1, 3)
    colors, counts = np.unique(flat, axis=0, return_counts=True)
    bg = colors[counts.argmax()]
    diff = np.abs(flat - bg).sum(axis=1)
    ink = (diff > 24).mean()  # >24 total-channel deviation = ink
    return float(ink)


def main() -> int:
    report = {"empty": [], "tiny": [], "missing": [], "ok": [], "unreadable": []}

    # 2026-09-06 audit fix: the previous dict built phantom keys
    # (f"{name}.png" pointing at main-dir paths that don't exist for the
    # numbered figures) which the exists() check silently skipped — the
    # glob-matched entries did the real work. Rebuilt below as a single
    # explicit list of (label, path) pairs, and the stray scan now covers
    # the supplementary directory too.
    expected: list[tuple[str, Path]] = []
    # numbered figures have varying suffixes; match by prefix
    for name in EXPECTED_MAIN:
        matches = sorted(FIG.glob(f"{name}_*.png"))
        if not matches:
            report["missing"].append(f"{name}_*.png (main)")
        else:
            expected.append((f"{name}", matches[0]))

    for name in EXPECTED_SUPP:
        p = SUP / f"{name}.png"
        if not p.exists():
            report["missing"].append(f"{name}.png (supplementary)")
        else:
            expected.append((name, p))

    audited = 0
    for key, path in expected:
        audited += 1
        size = path.stat().st_size
        name = Path(path).name
        if size < MIN_BYTES:
            report["tiny"].append({"file": name, "bytes": size})
            continue
        cov = _ink_coverage(path)
        if cov is None:
            report["unreadable"].append(name)
            continue
        if cov < INK_COVERAGE_MIN:
            report["empty"].append({"file": name, "ink_coverage": round(cov, 5)})
        else:
            report["ok"].append({"file": name, "ink_coverage": round(cov, 4)})

    # scan for stray non-numbered leftovers in BOTH main and supplementary
    stray_main = sorted(
        p.name for p in FIG.glob("*.png")
        if not any(p.name.startswith(f"{i:02d}_") for i in range(1, 34))
    )
    stray_supp = sorted(
        p.name for p in SUP.glob("*.png")
        if not p.name.startswith("fig_s")
    )
    report["stray_main_pngs"] = stray_main
    report["stray_supplementary_pngs"] = stray_supp
    # note: ink/byte thresholds CANNOT detect wrong-data or
    # truncation-biased content (e.g. the 2026-09-06 sensitivity-boxplot
    # bug passed as "ok"); content correctness is owned by the generating
    # scripts' own asserts, not by this file audit.
    report["scope_note"] = (
        "This audit verifies presence/size/ink only; semantic correctness "
        "of plotted data must be guaranteed by the generating scripts."
    )

    out_path = FIG / "figure_completeness_report.json"
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2)

    print(f"Audited {audited} figures (+{len(report['stray_main_pngs'])} stray).")
    if report["missing"]:
        print(f"  MISSING ({len(report['missing'])}):")
        for m in report["missing"]:
            print(f"    - {m}")
    if report["empty"]:
        print(f"  EMPTY  ({len(report['empty'])}):")
        for e in report["empty"]:
            print(f"    - {e['file']} (ink={e['ink_coverage']})")
    if report["tiny"]:
        print(f"  TINY   ({len(report['tiny'])}):")
        for t in report["tiny"]:
            print(f"    - {t['file']} ({t['bytes']} B)")
    if not (report["missing"] or report["empty"] or report["tiny"]):
        print("  All figures present and non-empty.")
    print(f"Report: {out_path}")

    return 0 if not (report["missing"] or report["empty"]) else 1


if __name__ == "__main__":
    sys.exit(main())
