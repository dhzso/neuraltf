#!/usr/bin/env python
"""Convert Fincher 2018 DGE to a subsampled h5ad (memory-efficient).

The Fincher atlas (GEO GSE111764) is provided as a single tab-separated
`PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz` file. The file
contains one row per cell, with the cell name in the first column followed
by integer counts for each gene (genes are listed in the header row, one
per column). This script subsamples cells to 10,000 (random seed 42) and
writes `datasets/processed/fincher_subsample.h5ad` so the NeuralTF pipeline
can load it without holding the full 26.5K-cell matrix in memory.

Usage:
    python scripts/convert_fincher.py [--cells 10000] [--seed 42]

If the raw DGE is not present, prints a clear message about how to download
it from GEO (accession GSE111764).
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
from scipy import sparse

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW = REPO_ROOT / "datasets" / "raw"
DGE = (
    RAW
    / "GSE111764_GEO_Fincher_atlas"
    / "GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz"
)
OUT_DIR = REPO_ROOT / "datasets" / "processed"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--cells",
        type=int,
        default=10000,
        help="Number of cells to subsample (default 10000; 0 = keep all)",
    )
    p.add_argument("--seed", type=int, default=42, help="RNG seed")
    p.add_argument(
        "--out",
        type=Path,
        default=OUT_DIR / "fincher_subsample.h5ad",
        help="Output h5ad path",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not DGE.exists():
        raise SystemExit(
            f"Fincher DGE not found at:\n  {DGE}\n\n"
            "Download from GEO (accession GSE111764):\n"
            "  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111764\n"
            f"and extract to:\n  {DGE.parent}"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)

    print("Reading Fincher DGE header...")
    with gzip.open(DGE, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    gene_ids = header[1:]
    n_genes = len(gene_ids)
    print(f"  Genes: {n_genes}")

    print("Counting cells...")
    cell_names: list[str] = []
    with gzip.open(DGE, "rt") as f:
        f.readline()
        for line in f:
            cell_names.append(line.split("\t", 1)[0].strip('"'))
    n_cells = len(cell_names)
    print(f"  Cells: {n_cells}")

    target = args.cells if args.cells and args.cells > 0 else n_cells
    if n_cells > target:
        rng = np.random.default_rng(args.seed)
        keep_idx = np.sort(rng.choice(n_cells, target, replace=False))
        keep_set = set(int(i) for i in keep_idx)
        print(f"  Subsampling to {target} cells")
    else:
        keep_idx = np.arange(n_cells)
        keep_set = set(int(i) for i in keep_idx)

    print("Building sparse matrix...")
    rows: list[int] = []
    cols: list[int] = []
    vals: list[int] = []
    kept_cell_names: list[str] = []

    with gzip.open(DGE, "rt") as f:
        f.readline()
        for i, line in enumerate(f):
            if i not in keep_set:
                continue
            parts = line.rstrip("\n").split("\t")
            cell_name = parts[0].strip('"')
            kept_cell_names.append(cell_name)
            row_idx = len(kept_cell_names) - 1
            counts = np.fromstring(
                "\t".join(parts[1:]), dtype=np.int32, sep="\t"
            )
            # Defend against extra columns. Many DGE files have a trailing
            # tab, producing one extra empty trailing element. We accept an
            # over-count of <= 1 (trailing-tab artefact: trim it) and pad
            # an under-count; any larger mismatch raises loudly so the
            # researcher knows the raw file orientation changed.
            if counts.size > n_genes:
                if counts.size == n_genes + 1:
                    # Trailing-tab artefact: trim the extra slot (np.fromstring
                    # would have produced 0 for the trailing empty string).
                    counts = counts[:n_genes]
                else:
                    raise ValueError(
                        f"Row {i} ({cell_name}) has {counts.size} count columns "
                        f"but the header declared {n_genes} genes. The raw file "
                        f"orientation may have changed."
                    )
            elif counts.size < n_genes:
                # Some cells may have fewer counts (pad with zeros).
                counts = np.pad(counts, (0, n_genes - counts.size))
            nz = np.flatnonzero(counts)
            if nz.size:
                rows.extend([row_idx] * nz.size)
                cols.extend(nz.tolist())
                vals.extend(counts[nz].tolist())

    print(f"  Kept cells: {len(kept_cell_names)}")
    print(f"  Non-zero entries: {len(vals)}")

    X = sparse.csr_matrix(
        (vals, (rows, cols)),
        shape=(len(kept_cell_names), n_genes),
        dtype=np.int32,
    )

    import anndata as ad

    adata = ad.AnnData(X=X)
    adata.var_names = gene_ids
    adata.obs_names = kept_cell_names
    adata.var["gene_id"] = gene_ids
    adata.uns["source"] = "GSE111764 (Fincher 2018)"
    adata.uns["subsample_cells"] = target
    adata.uns["subsample_seed"] = args.seed

    adata.write_h5ad(args.out)
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
