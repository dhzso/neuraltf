#!/usr/bin/env python
"""Convert Fincher 2018 DGE to a subsampled h5ad (memory-efficient).

The Fincher atlas (GEO GSE111764) is provided as a single tab-separated
`PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz` file. The first line
is the list of cell barcodes (one column per cell); every following line is
a gene ID followed by tab-separated counts (one row per gene, columns are
cells). By default the complete atlas is kept; pass --cells 10000 (seed 42)
to subsample for fast dev runs. Writes
`datasets/processed/fincher_subsample.h5ad` (cells x genes) so the NeuralTF
pipeline can load it without holding the full 50K-cell matrix in memory.

Usage:
    python scripts/convert_fincher.py [--cells 0] [--seed 42]

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
        default=0,
        help="Number of cells to keep (default 0 = the complete atlas; "
             "e.g. 10000 for a fast development run)",
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
    cell_names = [c.strip('"') for c in header]
    n_cells = len(cell_names)
    print(f"  Cells: {n_cells}")

    print("Counting genes...")
    gene_ids: list[str] = []
    with gzip.open(DGE, "rt") as f:
        f.readline()
        for line in f:
            gene_ids.append(line.split("\t", 1)[0].strip('"'))
    n_genes = len(gene_ids)
    print(f"  Genes: {n_genes}")

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
    kept_cell_names = [cell_names[i] for i in keep_idx]
    cell_pos = {int(c): i for i, c in enumerate(keep_idx)}

    with gzip.open(DGE, "rt") as f:
        f.readline()
        for gene_idx, line in enumerate(f):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            counts = parts[1:]
            if len(counts) > n_cells:
                # Trailing-tab artefact: trim the extra slot.
                if len(counts) == n_cells + 1:
                    counts = counts[:n_cells]
                else:
                    raise ValueError(
                        f"Gene row {gene_idx} ({gene_ids[gene_idx]}) has "
                        f"{len(counts)} count columns but the header declared "
                        f"{n_cells} cells. The raw file orientation may have "
                        f"changed."
                    )
            elif len(counts) < n_cells:
                counts = counts + [0] * (n_cells - len(counts))
            for cell_idx, tok in enumerate(counts):
                if cell_idx not in keep_set:
                    continue
                v = int(tok)
                if v > 0:
                    rows.append(gene_idx)
                    cols.append(cell_pos[cell_idx])
                    vals.append(v)

    print(f"  Kept cells: {len(kept_cell_names)}")
    print(f"  Non-zero entries: {len(vals)}")

    X = sparse.csr_matrix(
        (vals, (rows, cols)),
        shape=(n_genes, len(kept_cell_names)),
        dtype=np.int32,
    ).T
    print(f"  Matrix shape: {X.shape} (cells x genes)")

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
