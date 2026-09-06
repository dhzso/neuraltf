"""Tests for bioforge.omics.trajectory (PAGA + CellRank)."""
import warnings

import anndata as ad
import numpy as np
import pytest

from bioforge.omics.cluster import neighbors, pca, leiden
from bioforge.omics.normalize import normalize_total_log1p
from bioforge.omics.trajectory import paga, cellrank_terminal_states

_cellrank_available = True
try:
    import cellrank  # noqa: F401
except ImportError:
    _cellrank_available = False

needs_cellrank = pytest.mark.skipif(
    not _cellrank_available,
    reason="cellrank not installed (optional dep — requires native build)",
)

_scvelo_available = True
try:
    import scvelo  # noqa: F401
except ImportError:
    _scvelo_available = False

needs_scvelo = pytest.mark.skipif(
    not _scvelo_available,
    reason="scvelo not installed (optional dep — requires native build)",
)


def _prep_clustered(n_cells=200, n_genes=40, seed=0):
    rng = np.random.default_rng(seed)
    # Add some structure: two main groups with different means
    n_half = n_cells // 2
    x = np.zeros((n_cells, n_genes), dtype=np.float32)
    x[:n_half] = rng.poisson(1.5, (n_half, n_genes))
    x[n_half:] = rng.poisson(3.0, (n_cells - n_half, n_genes))
    adata = ad.AnnData(x.astype(np.float32))
    adata.obs_names = [f"c{i}" for i in range(n_cells)]
    adata.var_names = [f"g{j}" for j in range(n_genes)]
    normalize_total_log1p(adata)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pca(adata, n_comps=10, use_hvg=False)
        neighbors(adata, n_neighbors=10, n_pcs=10)
        leiden(adata, resolution=1.0)  # higher resolution for more clusters
    return adata


def test_paga_writes_connectivities() -> None:
    adata = _prep_clustered()
    paga(adata, groups="leiden")
    assert "paga" in adata.uns
    # PAGA returns connectivities and connectivities_tree
    assert "connectivities" in adata.uns["paga"]


def test_paga_unknown_groups_raises() -> None:
    adata = _prep_clustered()
    with pytest.raises(KeyError):
        paga(adata, groups="not_a_cluster_column")


@needs_scvelo
def test_velocity_without_spliced_unspliced_raises() -> None:
    from bioforge.omics.trajectory import velocity
    adata = _prep_clustered()
    with pytest.raises(KeyError, match="spliced"):
        velocity(adata)


@needs_cellrank
def test_cellrank_terminal_states_runs() -> None:
    """CellRank GPCCA needs an embedding + cluster_key column"""
    adata = _prep_clustered()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cellrank_terminal_states(adata, n_states=2, cluster_key="leiden")
    assert "terminal_states" in adata.obs
