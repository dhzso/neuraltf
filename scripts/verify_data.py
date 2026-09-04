#!/usr/bin/env python
"""Verify data integrity checksums for NeuralTF pipeline.

Usage:
    python scripts/verify_data.py              # verify all checksums
    python scripts/verify_data.py --generate   # generate checksums for current files
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RAW_DIR = REPO / "datasets" / "raw"

# Expected files with their relative paths â€” all pipeline-critical raw
# inputs (King mmc2-8, Perez MOESM5/19/22, Cui scRNA h5ad, Rosetta,
# go.obo included; MOESM2 is absent from the Perez download â€” see
# datasets/MANIFEST.md).
EXPECTED_FILES = [
    "GSE103633_GEO_Plass_atlas/GSE103633_RAW.tar",
    "GSE111764_GEO_Fincher_atlas/GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz",
    "OMIX003867_OMIX_Cui_atlas/OMIX003867-01/singlecell_h5ad/adata_scRNA_Annotated.h5ad",
    "Supplementary_Data_ King_2024/1-s2.0-S2211124724001712-mmc2.xlsx",
    "Supplementary_Data_ King_2024/1-s2.0-S2211124724001712-mmc3.xlsx",
    "Supplementary_Data_ King_2024/1-s2.0-S2211124724001712-mmc4.xlsx",
    "Supplementary_Data_ King_2024/1-s2.0-S2211124724001712-mmc5.xlsx",
    "Supplementary_Data_ King_2024/1-s2.0-S2211124724001712-mmc6.xlsx",
    "Supplementary_Data_ King_2024/1-s2.0-S2211124724001712-mmc7.xlsx",
    "Supplementary_Data_ Perez_2025/41467_2025_65712_MOESM5_ESM.xlsx",
    "Supplementary_Data_ Perez_2025/41467_2025_65712_MOESM19_ESM.xlsx",
    "Supplementary_Data_ Perez_2025/41467_2025_65712_MOESM22_ESM.xlsx",
    "smed_20140614.mapping.rosettastone.2020.txt",
    "go.obo",
]


def compute_sha256(path: Path) -> str:
    """Compute SHA256 hash of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description="Verify data checksums")
    parser.add_argument("--generate", action="store_true",
                        help="Generate checksums for current files (output JSON)")
    parser.add_argument("--manifest", action="store_true",
                        help="Print markdown table for MANIFEST.md")
    args = parser.parse_args()

    if args.generate:
        # Generate checksums
        checksums = {}
        for rel_path in EXPECTED_FILES:
            full_path = RAW_DIR / rel_path
            if full_path.exists():
                sha = compute_sha256(full_path)
                size = full_path.stat().st_size
                checksums[rel_path] = {"sha256": sha, "size_bytes": size}
                print(f"  {rel_path}: {sha} ({size:,} bytes)")
            else:
                checksums[rel_path] = {"sha256": "MISSING", "size_bytes": 0}
                print(f"  {rel_path}: MISSING")

        out_path = REPO / "datasets" / "checksums.json"
        out_path.write_text(json.dumps(checksums, indent=2))
        print(f"\nChecksums written to {out_path}")
        return 0

    if args.manifest:
        # Print markdown table for MANIFEST.md
        print("| File | SHA256 | Size | Source |")
        print("|------|--------|------|--------|")
        for rel_path in EXPECTED_FILES:
            full_path = RAW_DIR / rel_path
            if full_path.exists():
                sha = compute_sha256(full_path)
                size = full_path.stat().st_size
                size_str = f"{size/1e6:.1f} MB" if size > 1e6 else f"{size/1e3:.1f} KB"
                print(f"| `{rel_path}` | `{sha}` | {size_str} | GEO/Journal |")
            else:
                print(f"| `{rel_path}` | `MISSING` | â€” | â€” |")
        return 0

    # Default: verify against stored checksums
    checksums_path = REPO / "datasets" / "checksums.json"
    if not checksums_path.exists():
        print("No checksums.json found. Run with --generate first.")
        return 1

    expected = json.loads(checksums_path.read_text())
    all_ok = True

    print("=== Verifying Data Integrity ===")
    for rel_path, exp in expected.items():
        full_path = RAW_DIR / rel_path
        if not full_path.exists():
            print(f"  [FAIL] MISSING: {rel_path}")
            all_ok = False
            continue

        sha = compute_sha256(full_path)
        if sha == exp["sha256"]:
            print(f"  [PASS] {rel_path}")
        else:
            print(f"  [FAIL] MISMATCH: {rel_path}")
            print(f"    Expected: {exp['sha256']}")
            print(f"    Got:      {sha}")
            all_ok = False

    if all_ok:
        print("\n[PASS] All checksums verified.")
        return 0
    else:
        print("\n[FAIL] Some checksums failed.")
        return 1


if __name__ == "__main__":
    sys.exit(main())