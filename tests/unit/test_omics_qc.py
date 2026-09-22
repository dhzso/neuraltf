"""Tests for bioforge.omics.qc."""
import warnings

import anndata as ad
import numpy as np

from bioforge.omics.qc import (
    compute_qc_metrics,
    filter_cells,
    filter_genes,
    run_qc,
)


def _make_counts(n_cells: int = 500, n_genes: int = 300, mt: int = 5, seed: int = 0) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    # Use higher mean so min_counts=500 doesn't filter everything
    counts = rng.poisson(20, (n_cells, n_genes)).astype(np.float32)
    var_names = [f"g{i}" for i in range(n_genes)]
    if mt > 0:
        var_names[:mt] = [f"MT-{i}" for i in range(mt)]
    adata = ad.AnnData(counts)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = var_names
    return adata


def test_compute_qc_metrics_no_mt_prefix() -> None:
    adata = _make_counts()
    compute_qc_metrics(adata, mt_prefix=None)
    assert "total_counts" in adata.obs  # scanpy 1.11+ uses total_counts
    assert "n_genes_by_counts" in adata.obs


def test_compute_qc_metrics_with_mt_prefix() -> None:
    adata = _make_counts(mt=5)
    compute_qc_metrics(adata, mt_prefix="MT-")
    assert "mt" in adata.var
    assert adata.var["mt"].sum() == 5
    assert "pct_counts_mt" in adata.obs


def test_filter_cells_min_counts() -> None:
    adata = _make_counts()
    # Force the first 50 cells to have very few counts
    adata.X[:50] = 0
    filter_cells(adata, min_counts=10, min_genes=None, max_pct_mt=None)
    assert adata.n_obs < 500
    # Cells with all-zero counts are filtered. scanpy's filter_cells adds
    # 'n_counts' when compute_qc_metrics hasn't been called yet, or
    # 'total_counts' when it has. We accept either.
    counts_col = "n_counts" if "n_counts" in adata.obs else "total_counts"
    assert (adata.obs[counts_col] >= 10).all()


def test_filter_genes_min_cells() -> None:
    adata = _make_counts()
    # Force first 10 genes into 0 cells
    adata.X[:, :10] = 0
    n_before = adata.n_vars
    filter_genes(adata, min_cells=3)
    assert adata.n_vars < n_before
    # All remaining genes must have at least 3 cells with non-zero counts
    assert (np.asarray((adata.X > 0).sum(axis=0)).ravel() >= 3).all()


def test_run_qc_default_end_to_end() -> None:
    adata = _make_counts(mt=5, seed=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run_qc(adata, mt_prefix="MT-", max_pct_mt=20.0)
    assert adata.n_obs > 0
    assert adata.n_vars > 0
    assert "total_counts" in adata.obs


def test_run_qc_skip_mt_for_planarian() -> None:
    """Planarian transcripts have no 'MT-' prefix; QC must still run."""
    adata = _make_counts(mt=5, seed=2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run_qc(adata, mt_prefix=None, max_pct_mt=None)
    assert "pct_counts_mt" not in adata.obs
    assert "total_counts" in adata.obs
