"""Reader for Fincher et al. 2018 atlas files.

The Fincher atlas ships DGE matrices in
``datasets/raw/GSE111764_GEO_Fincher_atlas`` grouped by body region
(principal, brain, sexual clusterings). Gene IDs use **dd_Smed_v4**. This
module loads a chosen DGE file into an :class:`anndata.AnnData` so the rest
of the framework can plug it into the pipeline.

The reader only returns the matrix and its gene/cell names; cluster and
cell-type labels are not read here.
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from bioforge.core.logging import get_logger

logger = get_logger("evidence.readers.fincher")


def read_dge(path: str | Path) -> ad.AnnData:
    """Read a Fincher gzipped DGE matrix into an AnnData.

    The DGE txt files have gene IDs as rows and cell barcodes as columns.
    """
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    genes = df.index.astype(str)
    cells = df.columns.astype(str)
    X = df.to_numpy().T.astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.obs_names = cells
    adata.var_names = genes
    logger.info("loaded Fincher DGE: %d cells x %d genes from %s",
                adata.n_obs, adata.n_vars, Path(path).name)
    return adata
