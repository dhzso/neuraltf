#!/usr/bin/env python
"""Consolidate Plass 2018 per-cell DGE files into one subsampled h5ad.

The Plass atlas (GEO GSE103633) is provided as a tar archive of one
`<cell_name>.dge.txt.gz` file per cell (rows = gene ids, single count
column). This script streams the tar, builds a single cells x genes
sparse matrix, subsamples to 10,000 cells (random seed 42) and writes
`datasets/processed/plass_v6.h5ad`.

Usage:
    python scripts/consolidate_plass.py [--cells 10000] [--seed 42]

If the raw tar is not present, prints a clear message about how to download
it from GEO (accession GSE103633).
"""
from __future__ import annotations

import argparse
import gzip
import tarfile
from pathlib import Path

import numpy as np
from scipy import sparse

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW = REPO_ROOT / "datasets" / "raw"
TAR = RAW / "GSE103633_GEO_Plass_atlas" / "RAW.tar"
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
        default=OUT_DIR / "plass_v6.h5ad",
        help="Output h5ad path",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not TAR.exists():
        raise SystemExit(
            f"Plass RAW.tar not found at:\n  {TAR}\n\n"
            "Download from GEO (accession GSE103633):\n"
            "  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103633\n"
            f"and place RAW.tar at:\n  {TAR.parent}"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)

    print("Scanning Plass tar for cell DGE files...")
    with tarfile.open(TAR, "r") as tar:
        members = [m for m in tar.getmembers() if m.name.endswith(".dge.txt.gz")]
        print(f"  Found {len(members)} cell files")

        gene_to_idx: dict[str, int] = {}
        cell_names: list[str] = []
        rows: list[int] = []
        cols: list[int] = []
        vals: list[int] = []

        for i, member in enumerate(members):
            cell_name = member.name.replace(".dge.txt.gz", "").split("/")[-1]
            cell_names.append(cell_name)
            with tar.extractfile(member) as f:
                with gzip.open(f, "rt") as gz:
                    for line in gz:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 2:
                            continue
                        gene = parts[0]
                        try:
                            cnt = int(parts[1])
                        except ValueError:
                            continue
                        if cnt <= 0:
                            continue
                        if gene not in gene_to_idx:
                            gene_to_idx[gene] = len(gene_to_idx)
                        rows.append(i)
                        cols.append(gene_to_idx[gene])
                        vals.append(cnt)
            if (i + 1) % 2000 == 0:
                print(f"  Parsed {i + 1}/{len(members)} cells")

        gene_list = list(gene_to_idx.keys())
        print(f"  Total genes: {len(gene_list)}")

        matrix = sparse.csr_matrix(
            (vals, (rows, cols)),
            shape=(len(cell_names), len(gene_list)),
            dtype=np.int32,
        )
        print(f"  Matrix shape: {matrix.shape}, nnz: {matrix.nnz}")

    import anndata as ad

    adata = ad.AnnData(X=matrix)
    adata.var_names = gene_list
    adata.obs_names = cell_names
    adata.var["gene_id"] = gene_list
    adata.uns["source"] = "GSE103633 (Plass 2018)"

    target = args.cells if args.cells and args.cells > 0 else adata.n_obs
    if adata.n_obs > target:
        import scanpy as sc

        sc.pp.subsample(adata, n_obs=target, random_state=args.seed)
        adata.uns["subsample_cells"] = target
        adata.uns["subsample_seed"] = args.seed
        print(f"Subsampled to {adata.n_obs} cells")

    adata.write_h5ad(args.out)
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
