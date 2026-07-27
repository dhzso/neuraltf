"""Tests for bioforge.omics.batch (Harmony wrapper)."""
import warnings

import anndata as ad
import numpy as np
import pytest

from bioforge.omics.batch import run_harmony
from bioforge.omics.cluster import neighbors, pca
from bioforge.omics.normalize import normalize_total_log1p


def _prep_batches(n_per_batch=120, n_genes=50, seed=0):
    rng = np.random.default_rng(seed)
    a = rng.poisson(2, (2 * n_per_batch, n_genes)).astype(np.float32)
    adata = ad.AnnData(a)
    adata.obs_names = [f"c{i}" for i in range(2 * n_per_batch)]
    adata.var_names = [f"g{j}" for j in range(n_genes)]
    adata.obs["batch"] = ["a"] * n_per_batch + ["b"] * n_per_batch
    normalize_total_log1p(adata)
    return adata


def test_run_harmony_writes_corrected_obsm() -> None:
    adata = _prep_batches()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pca(adata, n_comps=10, use_hvg=False)
        run_harmony(adata, batch_key="batch", max_iter_harmony=3)
    assert "X_pca_harmony" in adata.obsm
    # Shape must match cells (no transpose bug)
    assert adata.obsm["X_pca_harmony"].shape == (adata.n_obs, 10)


def test_run_harmony_missing_basis_raises() -> None:
    adata = _prep_batches()
    # No PCA computed yet → run_harmony should raise KeyError
    with pytest.raises(KeyError, match="basis"):
        run_harmony(adata, batch_key="batch", basis="X_pca")


def test_run_harmony_missing_batch_key_raises() -> None:
    adata = _prep_batches()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pca(adata, n_comps=10, use_hvg=False)
    with pytest.raises(KeyError, match="batch key"):
        run_harmony(adata, batch_key="not_a_column")


def test_run_harmony_custom_basis_name() -> None:
    adata = _prep_batches()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pca(adata, n_comps=10, use_hvg=False)
        run_harmony(adata, batch_key="batch", max_iter_harmony=2,
                    corrected_basis="X_harmony_custom")
    assert "X_harmony_custom" in adata.obsm


def test_harmony_corrected_is_finite() -> None:
    adata = _prep_batches()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pca(adata, n_comps=10, use_hvg=False)
        run_harmony(adata, batch_key="batch", max_iter_harmony=3)
    assert np.isfinite(adata.obsm["X_pca_harmony"]).all()
