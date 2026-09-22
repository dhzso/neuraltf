"""Normalization and feature selection.

Mutates the AnnData in place and returns it. Always preserves raw counts
in ``adata.layers['counts']`` if they haven't already been stored, so
downstream callers can recover them.
"""

from __future__ import annotations

from typing import Optional

import scanpy as sc
from anndata import AnnData

from bioforge.core.logging import get_logger

logger = get_logger("omics.normalize")


def preserve_counts(adata: AnnData) -> AnnData:
    """Stash a copy of ``adata.X`` in ``layers['counts']`` if not already present.

    Idempotent — re-call on the same AnnData leaves existing counts intact.
    """
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
        logger.info("stashed raw counts in layers['counts']")
    return adata


def normalize_total_log1p(
    adata: AnnData,
    target_sum: float = 1e4,
) -> AnnData:
    """Apply scanpy ``normalize_total`` + ``log1p`` to ``adata.X``.

    Assumes the X matrix is the raw count matrix (it will preserve counts).
    """
    preserve_counts(adata)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    logger.info("normalized + log1p applied (target_sum=%g)", target_sum)
    return adata


def highly_variable_genes(
    adata: AnnData,
    n_top_genes: int = 2000,
    flavor: str = "seurat",
    batch_key: Optional[str] = None,
) -> AnnData:
    """Annotate highly variable genes (HVG). Mutates ``var['highly_variable']``."""
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor=flavor,
        batch_key=batch_key,
        subset=False,
    )
    n_hvg = int(adata.var["highly_variable"].sum())
    logger.info("flagged %d HVG (n_top=%d flavor=%s batch=%s)",
                n_hvg, n_top_genes, flavor, batch_key)
    return adata


def run_normalize(
    adata: AnnData,
    *,
    target_sum: float = 1e4,
    n_top_genes: int = 2000,
    hvg_flavor: str = "seurat",
    hvg_batch_key: Optional[str] = None,
) -> AnnData:
    """Convenience: normalize + HVG in one step."""
    normalize_total_log1p(adata, target_sum=target_sum)
    highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor=hvg_flavor,
        batch_key=hvg_batch_key,
    )
    logger.info("run_normalize complete")
    return adata


__all__ = [
    "preserve_counts",
    "normalize_total_log1p",
    "highly_variable_genes",
    "run_normalize",
]
