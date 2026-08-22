"""Quality control for single-cell count matrices.

Provides opinionated QC functions that mutate the AnnData in place and
return the same object. Designed for droplet-based scRNA-seq (Drop-seq,
10x). Mitochondrial percent is computed when ``mt_prefix`` is supplied;
otherwise the QC metric is skipped (the planarian transcriptome has no
mitochondrial chromosome, so planarian analyses may skip it).
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import scanpy as sc
from anndata import AnnData

from bioforge.core.logging import get_logger

logger = get_logger("omics.qc")


def compute_qc_metrics(
    adata: AnnData,
    mt_prefix: Optional[str] = "MT-",
) -> AnnData:
    """Compute standard QC metrics, including mitochondrial percent.

    Parameters
    ----------
    adata
        AnnData with raw counts in ``adata.X`` (or ``.layers['counts']``).
    mt_prefix
        Prefix for mitochondrial genes (e.g. ``"MT-"`` in human GRCh38).
        Set to ``None`` to skip the mitochondrial flag (planarian use).

    Returns
    -------
    AnnData
        The same object with ``obs['total_counts']``, ``obs['n_genes_by_counts']``
        (and ``obs['pct_counts_mt']`` when mt_prefix is given) populated.
    """
    logger.info(
        "computing qc metrics on %d cells x %d genes (mt_prefix=%s)",
        adata.n_obs, adata.n_vars, mt_prefix,
    )
    if mt_prefix:
        _flag_mt_genes(adata, mt_prefix)
        sc.pp.calculate_qc_metrics(
            adata,
            percent_top=None,
            log1p=False,
            inplace=True,
            qc_vars=["mt"],
        )
    else:
        sc.pp.calculate_qc_metrics(
            adata,
            percent_top=None,
            log1p=False,
            inplace=True,
        )
    return adata


def _flag_mt_genes(adata: AnnData, prefix: str) -> bool:
    """Tag mitochondrial genes via a ``var['mt']`` boolean Series.

    Returns True if at least one gene matched the prefix; otherwise returns
    False so callers know the flag is empty.
    """
    adata.var["mt"] = adata.var_names.str.startswith(prefix)
    n_mt = int(adata.var["mt"].sum())
    logger.info("flagged %d mitochondrial genes with prefix '%s'", n_mt, prefix)
    return n_mt > 0


def filter_cells(
    adata: AnnData,
    min_counts: int = 500,
    min_genes: Optional[int] = 200,
    max_pct_mt: Optional[float] = 20.0,
) -> AnnData:
    """Filter low-quality cells by counts, gene count, and mitochondrial %.

    Mutates ``adata`` in place where possible (the total-counts filter uses
    scanpy's in-place routine; subsequent gene/mt-percent filters use
    boolean masking). Recomputes ``n_genes_by_counts`` reference metrics are
    left as scanpy set them.

    Parameters
    ----------
    adata
        Input AnnData.
    min_counts
        Minimum total UMI counts per cell.
    min_genes
        Minimum number of genes detected per cell. Set to ``None`` to skip.
    max_pct_mt
        Maximum mitochondrial percent allowed. Skipped if ``None`` or if
        ``pct_counts_mt`` is not in ``adata.obs``.
    """
    n_before = adata.n_obs

    # Total-count filter (scanpy, in place)
    sc.pp.filter_cells(adata, min_counts=min_counts)

    # Gene-count filter (apply mask and replace contents via public API)
    if min_genes is not None and "n_genes_by_counts" in adata.obs:
        mask = (adata.obs["n_genes_by_counts"] >= min_genes).to_numpy()
        if not mask.all():
            _inplace_subset_obs(adata, mask)

    # Mitochondrial percent filter
    if max_pct_mt is not None and "pct_counts_mt" in adata.obs:
        mask = (adata.obs["pct_counts_mt"] < max_pct_mt).to_numpy()
        if not mask.all():
            _inplace_subset_obs(adata, mask)

    logger.info(
        "filter_cells: kept %d/%d cells (min_counts=%d min_genes=%s max_pct_mt=%s)",
        adata.n_obs, n_before, min_counts, min_genes, max_pct_mt,
    )
    return adata


def _inplace_subset_obs(adata: AnnData, mask) -> None:
    """Subset ``adata`` observations in place using a boolean mask.

    Replaces the AnnData's underlying storage so callers passing the object
    by reference see the filtered state without an extra copy step.
    """
    new = adata[mask].copy()
    adata._init_as_actual(new)


def filter_genes(
    adata: AnnData,
    min_cells: int = 3,
) -> AnnData:
    """Filter genes expressed in fewer than ``min_cells`` cells."""
    n_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    logger.info("filter_genes: kept %d/%d genes (min_cells=%d)",
                adata.n_vars, n_before, min_cells)
    return adata


def run_qc(
    adata: AnnData,
    *,
    min_counts_per_cell: int = 500,
    min_genes_per_cell: int = 200,
    min_cells_per_gene: int = 3,
    mt_prefix: Optional[str] = "MT-",
    max_pct_mt: Optional[float] = 20.0,
) -> AnnData:
    """Run the standard BioForge QC pipeline (metrics + filters).

    Convenience entry point that chains :func:`compute_qc_metrics`,
    :func:`filter_cells`, and :func:`filter_genes`.
    """
    compute_qc_metrics(adata, mt_prefix=mt_prefix)
    filter_cells(
        adata,
        min_counts=min_counts_per_cell,
        min_genes=min_genes_per_cell,
        max_pct_mt=max_pct_mt,
    )
    filter_genes(adata, min_cells=min_cells_per_gene)
    logger.info("run_qc complete: %d cells x %d genes", adata.n_obs, adata.n_vars)
    return adata


__all__ = [
    "compute_qc_metrics",
    "filter_cells",
    "filter_genes",
    "run_qc",
]
