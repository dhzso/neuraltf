#!/usr/bin/env python
"""Consolidate Plass 2018 per-cell DGE files into one subsampled h5ad.

The Plass atlas (GEO GSE103633) is provided as a tar archive containing one
`<sample>_DGE_CLEAN.MULTI.txt.gz` file per sample (11 samples for this
study). Each such file is a cells x genes DGE table: the first line is
`GENE<tab><cell barcode>...`, and every following line is a gene id followed
by tab-separated counts (gene x cell). This script streams the tar, joins all
samples into a single cells x genes sparse matrix (cell barcodes are prefixed
with their sample tag to keep them unique across samples), subsamples to
10,000 cells (random seed 42) and writes `datasets/processed/plass_v6.h5ad`.

Usage:
    python scripts/consolidate_plass.py [--cells 10000] [--seed 42]
    python scripts/consolidate_plass.py [--tar datasets/raw/.../GSE103633_RAW.tar]

If a `--tar` path is not given, the tar is located automatically under
`datasets/raw` (it may be named RAW.tar or GSE103633_RAW.tar); otherwise a
clear message about how to download it from GEO (accession GSE103633) is
printed.
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
OUT_DIR = REPO_ROOT / "datasets" / "processed"

# GEO names the supplementary tar "GSE103633_RAW.tar". People usually drop it,
# plus the other GSE103633 files, into a single folder, so we do not rely on a
# single hardcoded location.
_PLASS_TAR_CANDIDATES = (
    RAW / "Plass_2018" / "RAW.tar",
    RAW / "GSE103633_GEO_Plass_atlas" / "GSE103633_RAW.tar",
    RAW / "GSE103633_GEO_Plass_atlas" / "RAW.tar",
    RAW / "GSE103633_RAW.tar",
    RAW / "Plass_2018" / "GSE103633_RAW.tar",
)

# Member files inside the tar that carry the DGE tables (case-insensitive).
_DGE_SUFFIXES = ("DGE_CLEAN.MULTI.TXT.GZ", "DGE.TXT.GZ")


def resolve_plass_tar() -> Path:
    """Locate the Plass tar wherever it was downloaded to.

    Tries the known locations first, then falls back to searching
    ``datasets/raw`` for any ``*.tar`` whose name contains "RAW" or
    "GSE103633". Raises SystemExit with a clear message if nothing matches.
    """
    for path in _PLASS_TAR_CANDIDATES:
        if path.is_file():
            return path

    hits: list[Path] = []
    for path in sorted(RAW.rglob("*")):
        if not path.is_file() or path.suffix.lower() != ".tar":
            continue
        base = path.name.upper()
        if "RAW" in base or "GSE103633" in base:
            hits.append(path)

    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        # Prefer the exact names GEO produces (GSE103633_RAW.tar) or the
        # plain "RAW.tar" the docs describe.
        for path in hits:
            if path.name.upper() in {"GSE103633_RAW.TAR", "RAW.TAR"}:
                return path
        raise SystemExit(
            "Multiple candidate Plass tar files found:\n"
            + "\n".join(f"  {p}" for p in hits)
            + "\nRename the correct one (the tar of per-sample DGE files) "
            "to RAW.tar or GSE103633_RAW.tar."
        )

    raise SystemExit(
        "Plass RAW.tar not found under datasets/raw/.\n"
        "Download from GEO (accession GSE103633):\n"
        "  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103633\n"
        "and place the downloaded tar (it will be named GSE103633_RAW.tar, "
        "the one containing the per-sample DGE .txt.gz files) somewhere under:\n"
        "  datasets/raw/"
    )


def _member_suffix(name: str) -> str | None:
    """Return the known DGE suffix that ``name`` ends with, else None."""
    upper = name.upper()
    for suffix in _DGE_SUFFIXES:
        if upper.endswith(suffix):
            return suffix
    return None


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
        "--tar",
        type=Path,
        default=None,
        help="Explicit path to the Plass tar (default: auto-locate)",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=OUT_DIR / "plass_v6.h5ad",
        help="Output h5ad path",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    tar_path = args.tar if args.tar is not None else resolve_plass_tar()
    print(f"Using Plass tar: {tar_path}")

    args.out.parent.mkdir(parents=True, exist_ok=True)

    all_cells: list[str] = []
    gene_to_idx: dict[str, int] = {}
    rows: list[int] = []
    cols: list[int] = []
    vals: list[int] = []

    with tarfile.open(tar_path, "r") as tar:
        dge_members = [m for m in tar.getmembers() if _member_suffix(m.name)]
        print(f"  Found {len(dge_members)} sample DGE files")
        if not dge_members:
            raise SystemExit(
                "No sample DGE files (names ending in "
                "DGE_CLEAN.MULTI.txt.gz or DGE.txt.gz) found inside the tar. "
                "Confirm you downloaded the RAW data tar from GSE103633."
            )

        for member in dge_members:
            suffix = _member_suffix(member.name)
            tag = member.name[:-len(suffix)].replace("\\", "/").split("/")[-1]
            if not tag:
                tag = member.name

            base_cell_idx = len(all_cells)
            with tar.extractfile(member) as f, gzip.open(f, "rt") as gz:
                header = gz.readline().rstrip("\n").split("\t")
                if len(header) < 2 or header[0].upper() != "GENE":
                    print(f"  !! {member.name}: unexpected header, skipping")
                    continue
                n_cells = len(header) - 1
                for bc in header[1:]:
                    all_cells.append(f"{tag}__{bc}")
                print(f"  {tag}: {n_cells} cells")

                for line in gz:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2 or parts[0].strip() == "":
                        continue
                    gene = parts[0]
                    counts = parts[1:]
                    if len(counts) > n_cells:
                        if len(counts) == n_cells + 1:
                            counts = counts[:n_cells]
                        else:
                            raise ValueError(
                                f"{member.name} row '{gene}' has {len(counts)} "
                                f"count columns but the header declared "
                                f"{n_cells} cells. Raw orientation changed?"
                            )
                    elif len(counts) < n_cells:
                        counts = counts + [0] * (n_cells - len(counts))

                    if gene not in gene_to_idx:
                        gene_to_idx[gene] = len(gene_to_idx)
                    gidx = gene_to_idx[gene]
                    for j, tok in enumerate(counts):
                        v = int(tok)
                        if v > 0:
                            rows.append(base_cell_idx + j)
                            cols.append(gidx)
                            vals.append(v)

    n_cells = len(all_cells)
    n_genes = len(gene_to_idx)
    print(f"  Total cells: {n_cells}, total genes: {n_genes}")
    if n_cells == 0 or n_genes == 0:
        raise SystemExit(
            "No expression data found in the tar. Confirm the tar contains "
            "the per-sample DGE_CLEAN.MULTI.txt.gz files from GSE103633."
        )

    gene_list = list(gene_to_idx)
    matrix = sparse.csr_matrix(
        (vals, (rows, cols)),
        shape=(n_cells, n_genes),
        dtype=np.int64,
    )
    print(f"  Matrix shape: {matrix.shape}, nnz: {matrix.nnz}")

    import anndata as ad

    adata = ad.AnnData(X=matrix)
    adata.var_names = gene_list
    adata.obs_names = all_cells
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