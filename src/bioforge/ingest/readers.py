"""Concrete matrix readers. Each returns an :class:`anndata.AnnData`."""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Union

import anndata as ad
import numpy as np
import pandas as pd

from bioforge.core.logging import get_logger

logger = get_logger("ingest.readers")


def read_dge_gz(path: Union[str, Path]) -> ad.AnnData:
    """Read a gzipped DGE-style TSV (genes as rows, cells as columns)."""
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    return _df_genes_to_anndata(df, path)


def read_dge_txt(path: Union[str, Path]) -> ad.AnnData:
    """Read an uncompressed DGE-style TSV."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    return _df_genes_to_anndata(df, path)


def read_h5ad(path: Union[str, Path]) -> ad.AnnData:
    return ad.read_h5ad(path)


def read_10x_mtx(path: Union[str, Path], **kwargs) -> ad.AnnData:
    """Read a 10x matrix.mtx + barcodes + features directory."""
    import scanpy as sc
    return sc.read_10x_mtx(str(path), **kwargs)


def read_csv_matrix(path: Union[str, Path]) -> ad.AnnData:
    """Read a plain CSV matrix (genes as rows by default)."""
    df = pd.read_csv(path, index_col=0)
    return _df_genes_to_anndata(df, path)


def read_tsv_matrix(path: Union[str, Path]) -> ad.AnnData:
    """Read a tab-separated matrix (genes as rows)."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    return _df_genes_to_anndata(df, path)


def _df_genes_to_anndata(df: pd.DataFrame, path: Union[str, Path]) -> ad.AnnData:
    genes = df.index.astype(str)
    cells = df.columns.astype(str)
    X = df.to_numpy().T.astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.obs_names = cells
    adata.var_names = genes
    adata.uns["_format_source_file"] = str(path)
    logger.info(
        "loaded %s: %d cells x %d genes",
        Path(path).name, adata.n_obs, adata.n_vars,
    )
    return adata


def peek_first_bytes(path: Union[str, Path], n: int = 512) -> bytes:
    """Read up to ``n`` bytes from a file or gzipped file (transparent)."""
    with open(path, "rb") as fh:
        head = fh.read(n)
    if head.startswith(b"\x1f\x8b"):  # gzip magic
        # Decompress just the leading chunk for sniffing.
        with gzip.open(path, "rb") as gf:
            head = gf.read(n)
    return head
