"""Tests for bioforge.omics.normalize."""
import anndata as ad
import numpy as np

from bioforge.omics.normalize import (
    highly_variable_genes,
    normalize_total_log1p,
    preserve_counts,
    run_normalize,
)


def _make_log_input(n_cells=300, n_genes=80, seed=0):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(2, (n_cells, n_genes)).astype(np.float32)
    adata = ad.AnnData(counts)
    adata.obs_names = [f"c{i}" for i in range(n_cells)]
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    return adata


def test_preserve_counts_idempotent() -> None:
    adata = _make_log_input()
    preserve_counts(adata)
    preserve_counts(adata)  # second call should not stomp the stored copy
    assert "counts" in adata.layers
    np.testing.assert_allclose(adata.layers["counts"], adata.X)


def test_normalize_total_log1p_changes_X() -> None:
    adata = _make_log_input()
    X_raw = adata.X.copy()
    normalize_total_log1p(adata, target_sum=1e4)
    # After log1p + normalize, X is no longer equal to raw counts
    assert not np.allclose(adata.X, X_raw)
    # And raw counts are preserved in layers['counts']
    np.testing.assert_allclose(adata.layers["counts"], X_raw)


def test_highly_variable_genes_flag() -> None:
    adata = _make_log_input()
    normalize_total_log1p(adata)
    highly_variable_genes(adata, n_top_genes=20, flavor="seurat")
    assert "highly_variable" in adata.var
    # n_top_genes is an upper bound in scanpy; allow ±2
    n_hvg = int(adata.var["highly_variable"].sum())
    assert 18 <= n_hvg <= 22


def test_run_normalize_pipeline() -> None:
    adata = _make_log_input()
    run_normalize(adata, n_top_genes=15)
    assert "counts" in adata.layers
    assert "highly_variable" in adata.var
    n_hvg = int(adata.var["highly_variable"].sum())
    assert 13 <= n_hvg <= 17
