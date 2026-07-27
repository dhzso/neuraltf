"""Trajectory inference: PAGA, RNA velocity, CellRank helpers.

These thin wrappers expose the standard trajectory tools used in planarian
scRNA-seq analysis (Plass et al. 2018, King et al. 2024) while shielding
BioForge callers from differences across scverse minor versions.
"""

from __future__ import annotations

from typing import Iterable, Optional

import scanpy as sc
import scvelo as scv
from anndata import AnnData

from bioforge.core.logging import get_logger

logger = get_logger("omics.trajectory")


def paga(
    adata: AnnData,
    groups: str = "leiden",
    use_rna_velocity: bool = False,
) -> AnnData:
    """Run PAGA (partition-based graph abstraction).

    Parameters
    ----------
    groups
        ``obs`` column with the cluster labels to abstract over
        (default ``'leiden'``).
    use_rna_velocity
        If True, use the RNA velocity-derived transition matrix if available
        (requires scVelo pre-computed). Otherwise use the standard kNN graph.
    """
    sc.tl.paga(adata, groups=groups, use_rna_velocity=use_rna_velocity)
    logger.info("PAGA computed over groups='%s' (use_rna_velocity=%s)",
                groups, use_rna_velocity)
    return adata


def velocity(
    adata: AnnData,
    mode: str = "stochastic",
    n_jobs: int = 1,
) -> AnnData:
    """Run scVelo RNA velocity.

    Expects spliced/unspliced layers to be present in ``adata.layers``
    (typically loaded via scv.read_loom or your own loader). Lookups for
    spliced/unspliced are exactly as scVelo expects them.
    """
    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        raise KeyError(
            "RNA velocity requires adata.layers['spliced'] and "
            "adata.layers['unspliced']. Load a loom file with scv.read_loom "
            "or fill these layers manually."
        )
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(adata, mode=mode, n_jobs=n_jobs)
    scv.tl.velocity_graph(adata)
    logger.info("RNA velocity computed (mode=%s, n_jobs=%d)", mode, n_jobs)
    return adata


def cellrank_terminal_states(
    adata: AnnData,
    n_states: int = 4,
    cluster_key: str = "leiden",
) -> AnnData:
    """Initialize CellRank and compute terminal macrostates using GPCCA.

    Mutates adata with CellRank outputs (``obs['terminal_states']``,
    ``obsm['macrostates']``, etc.). Returns the AnnData for chaining.
    """
    import cellrank as cr
    from cellrank.kernels import VelocityKernel, ConnectivityKernel
    cr.settings.verbosity = 0

    # Prefer velocity-based transition matrix if velocity is available
    if "velocity_graph" in adata.obsp or "velocyto_transitions" in adata.obsp:
        kernel = VelocityKernel(adata)
    else:
        kernel = ConnectivityKernel(adata)
    kernel.compute_transition_matrix()
    estimator = cr.estimators.GPCCA(kernel)
    estimator.compute_schur(n_components=min(n_states * 2, 20))
    estimator.compute_macrostates(n_states=n_states, cluster_key=cluster_key)
    estimator.set_terminal_states(n_cells=30, cluster_key=cluster_key)
    estimator.compute_fate_probabilities()
    # Write canonical CellRank outputs back into the AnnData for downstream code
    adata.obs["terminal_states"] = estimator.terminal_states
    if estimator.fate_probabilities is not None:
        adata.obsm["lineage_protocol"] = estimator.fate_probabilities.X
    logger.info("CellRank terminal states: %d (cluster_key=%s)", n_states, cluster_key)
    return adata


__all__ = [
    "paga",
    "velocity",
    "cellrank_terminal_states",
]
