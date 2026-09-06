"""Tests for bioforge.omics.cluster."""
import warnings

import anndata as ad
import numpy as np

from bioforge.omics.cluster import (
    leiden,
    neighbors,
    pca,
    run_cluster,
    umap,
)
from bioforge.omics.normalize import normalize_total_log1p


def _prep(n_cells=250, n_genes=60, seed=1):
    rng = np.random.default_rng(seed)
    adata = ad.AnnData(rng.poisson(2, (n_cells, n_genes)).astype(np.float32))
    adata.obs_names = [f"c{i}" for i in range(n_cells)]
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    normalize_total_log1p(adata)
    return adata


def test_pca_writes_obsm(tmp_path) -> None:
    adata = _prep()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pca(adata, n_comps=10, use_hvg=False)
    assert "X_pca" in adata.obsm
    assert adata.obsm["X_pca"].shape == (adata.n_obs, 10)


def test_neighbors_writes_obsp() -> None:
    adata = _prep()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pca(adata, n_comps=10, use_hvg=False)
        neighbors(adata, n_neighbors=10, n_pcs=10)
    assert "neighbors" in adata.uns
    assert "distances" in adata.obsp


def test_leiden_writes_obs_column() -> None:
    adata = _prep()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pca(adata, n_comps=10, use_hvg=False)
        neighbors(adata, n_neighbors=10, n_pcs=10)
        leiden(adata, resolution=0.5)
    assert "leiden" in adata.obs
    assert adata.obs["leiden"].nunique() >= 1


def test_umap_writes_obsm() -> None:
    adata = _prep()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pca(adata, n_comps=10, use_hvg=False)
        neighbors(adata, n_neighbors=10, n_pcs=10)
        umap(adata)
    assert "X_umap" in adata.obsm
    assert adata.obsm["X_umap"].shape == (adata.n_obs, 2)


def test_run_cluster_full_pipeline() -> None:
    adata = _prep()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run_cluster(adata, n_comps=15, n_neighbors=10, n_pcs=10,
                    resolution=0.5, use_hvg=False, compute_umap=True)
    assert "X_pca" in adata.obsm
    assert "leiden" in adata.obs
    assert "X_umap" in adata.obsm


def test_custom_leiden_key_added() -> None:
    adata = _prep()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pca(adata, n_comps=10, use_hvg=False)
        neighbors(adata, n_neighbors=10, n_pcs=10)
        leiden(adata, resolution=1.0, key_added="cluster_v1")
    assert "cluster_v1" in adata.obs
