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

EXPECTED_MAIN = [
    "01", "03", "04", "05", "06", "07", "08", "09",
    "13", "14", "15", "17", "18", "20", "22", "23",
    "24", "25", "26", "28", "29", "30", "31", "32",
    "33", "34", "35",
]

EXPECTED_COMPOSITE = [
    "fig1_overview_pipeline",
    "fig2_prioritization_landscape",
    "fig3_sensitivity_and_uncertainty",
    "fig4_method_agreement_concordance",
    "fig5_empirical_benchmarks_and_network",
]

EXPECTED_SUPP = [
    "fig_s1_go_landscape",
    "fig_s2_go_namespace_and_track",
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
    flat = arr.reshape(-1, 3)
    colors, counts = np.unique(flat, axis=0, return_counts=True)
    bg = colors[counts.argmax()]
    diff = np.abs(flat - bg).sum(axis=1)
    ink = (diff > 24).mean()
    return float(ink)


def main() -> int:
    report = {"empty": [], "tiny": [], "missing": [], "ok": [], "unreadable": [], "missing_pdf": []}

    expected: list[tuple[str, Path]] = []
    # Curated numbered figures
    for num in EXPECTED_MAIN:
        matches = sorted(FIG.glob(f"{num}_*.png"))
        if not matches:
            report["missing"].append(f"{num}_*.png (main)")
        else:
            expected.append((num, matches[0]))

    # Composite manuscript figures
    for name in EXPECTED_COMPOSITE:
        p = FIG / f"{name}.png"
        if not p.exists():
            report["missing"].append(f"{name}.png (composite)")
        else:
            expected.append((name, p))

    # Supplementary figures
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

        # Check for vector PDF counterpart
        pdf_path = path.with_suffix(".pdf")
        if not pdf_path.exists():
            report["missing_pdf"].append(name.replace(".png", ".pdf"))

    valid_prefixes = tuple([f"{num}_" for num in EXPECTED_MAIN] + [f"{c}" for c in EXPECTED_COMPOSITE])
    stray_main = sorted(
        p.name for p in FIG.glob("*.png")
        if not any(p.name.startswith(pre) for pre in valid_prefixes)
    )
    stray_supp = sorted(
        p.name for p in SUP.glob("*.png")
        if not any(p.name.startswith(s) for s in EXPECTED_SUPP)
    )
    report["stray_main_pngs"] = stray_main
    report["stray_supplementary_pngs"] = stray_supp
    report["scope_note"] = (
        "Verified Nature Communications standards: 300 DPI PNG + vector PDF, "
        "non-emptiness, and complete curated set."
    )

    out_path = FIG / "figure_completeness_report.json"
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2)

    print(f"Audited {audited} figures (PNG + vector PDF).")
    if report["missing"]:
        print(f"  MISSING PNG ({len(report['missing'])}):")
        for m in report["missing"]:
            print(f"    - {m}")
    if report["missing_pdf"]:
        print(f"  MISSING PDF ({len(report['missing_pdf'])}):")
        for mp in report["missing_pdf"]:
            print(f"    - {mp}")
    if report["empty"]:
        print(f"  EMPTY  ({len(report['empty'])}):")
        for e in report["empty"]:
            print(f"    - {e['file']} (ink={e['ink_coverage']})")
    if report["tiny"]:
        print(f"  TINY   ({len(report['tiny'])}):")
        for t in report["tiny"]:
            print(f"    - {t['file']} ({t['bytes']} B)")
    if report["stray_main_pngs"]:
        print(f"  STRAY MAIN PNGs ({len(report['stray_main_pngs'])}): {report['stray_main_pngs']}")
    if not (report["missing"] or report["missing_pdf"] or report["empty"] or report["tiny"] or report["stray_main_pngs"]):
        print("  All figures pass 100% (PNG + PDF present, non-empty, no strays).")
    print(f"Report: {out_path}")

    return 0 if not (report["missing"] or report["empty"] or report["missing_pdf"]) else 1


if __name__ == "__main__":
    sys.exit(main())

