"""Reader for Plass et al. 2018 atlas files.

The Plass atlas ships as ``GSE103633_RAW.tar`` (a tarball of per-cell
expression files) plus the series matrix and a contigs FASTA. Gene IDs
use **dd_Smed_v6**. The reader focus here mirrors :mod:`fincher`: a thin
loader that exposes cells × genes data to the framework.

For the thesis we expect to pre-extract one usable matrix file; the
reader supports either an h5ad path or a gzipped DGE-style txt.
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from bioforge.core.logging import get_logger

logger = get_logger("evidence.readers.plass")


def read_plass_matrix(path: str | Path) -> ad.AnnData:
    """Load a Plass atlas matrix file into an AnnData.

    Supports two file extensions:
    - ``.h5ad`` — read directly with :func:`anndata.read_h5ad`.
    - ``.txt.gz`` — TSV with genes as rows, cells as columns (Fincher-like).
    """
    p = Path(path)
    if p.suffix == ".h5ad":
        adata = ad.read_h5ad(p)
        logger.info("loaded Plass h5ad: %d cells x %d genes", adata.n_obs, adata.n_vars)
        return adata
    if "".join(p.suffixes[-2:]) == ".txt.gz":
        df = pd.read_csv(p, sep="\t", index_col=0, compression="gzip")
        genes = df.index.astype(str)
        cells = df.columns.astype(str)
        X = df.to_numpy().T.astype(np.float32)
        adata = ad.AnnData(X=X)
        adata.obs_names = cells
        adata.var_names = genes
        logger.info("loaded Plass DGE: %d cells x %d genes", adata.n_obs, adata.n_vars)
        return adata
    raise ValueError(f"unsupported Plass matrix extension: {p.name}")
