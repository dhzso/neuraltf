"""Clustering: PCA, neighborhood graph, Leiden, UMAP."""

from __future__ import annotations

from typing import Optional

import scanpy as sc
from anndata import AnnData

from bioforge.core.logging import get_logger

logger = get_logger("omics.cluster")


def pca(adata: AnnData, n_comps: int = 50, use_hvg: bool = True) -> AnnData:
    """Run PCA on the (normalized log) matrix.

    When ``use_hvg=True``, only ``var['highly_variable']==True`` genes are
    used; requires normalize module has flagged HVGs.
    """
    mask_var = "highly_variable" if (use_hvg and "highly_variable" in adata.var) else None
    sc.pp.pca(
        adata,
        n_comps=n_comps,
        mask_var=mask_var,
    )
    logger.info("PCA computed: n_comps=%d mask_var=%s", n_comps, mask_var)
    return adata



def neighbors(
    adata: AnnData,
    n_neighbors: int = 15,
    n_pcs: int = 30,
    use_rep: Optional[str] = None,
) -> AnnData:
    """Compute the kNN graph used by clustering and embeddings."""
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=use_rep,
    )
    logger.info(
        "neighbors graph built: n_neighbors=%d n_pcs=%d use_rep=%s",
        n_neighbors, n_pcs, use_rep,
    )
    return adata


def leiden(
    adata: AnnData,
    resolution: float = 1.0,
    flavor: str = "igraph",
    directed: bool = False,
    key_added: str = "leiden",
    use_rep: Optional[str] = None,
) -> AnnData:
    """Leiden community detection.

    Notes
    -----
    Default ``flavor='igraph'`` matches scanpy 1.11+'s recommended default
    to avoid the deprecation of ``flavor='leidenalg'`` behavior in future
    versions.
    """
    sc.tl.leiden(
        adata,
        resolution=resolution,
        flavor=flavor,
        directed=directed,
        key_added=key_added,
        neighbors_key=None,
    )
    n_cl = int(adata.obs[key_added].nunique())
    logger.info("leiden clusters: %d (resolution=%g)", n_cl, resolution)
    return adata


def umap(adata: AnnData) -> AnnData:
    """Run UMAP using the existing neighbor graph."""
    sc.tl.umap(adata)
    logger.info("UMAP embedding computed")
    return adata


def run_cluster(
    adata: AnnData,
    *,
    n_comps: int = 50,
    n_neighbors: int = 15,
    n_pcs: int = 30,
    resolution: float = 1.0,
    use_hvg: bool = True,
    compute_umap: bool = True,
) -> AnnData:
    """End-to-end clustering pipeline: PCA → neighbors → leiden → UMAP."""
    pca(adata, n_comps=n_comps, use_hvg=use_hvg)
    neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    leiden(adata, resolution=resolution)
    if compute_umap:
        umap(adata)
    logger.info("run_cluster complete")
    return adata


__all__ = ["pca", "neighbors", "leiden", "umap", "run_cluster"]
