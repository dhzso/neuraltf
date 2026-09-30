"""Batch correction via Harmony.

Wraps harmonypy directly instead of using scanpy's
``sc.external.pp.harmony_integrate``: on scanpy 1.11.5 that helper transposes
``ho.Z_corr``, which harmonypy 2.0 already returns as (cells × comps).

harmonypy is an *optional* dependency (it requires a native BLAS build
that is painful on a fresh Windows install). The module imports cleanly
without it; ``run_harmony()`` raises a clear ``ImportError`` if invoked
when harmonypy isn't installed.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
from anndata import AnnData

from bioforge.core.logging import get_logger

logger = get_logger("omics.batch")


def _import_harmonypy():
    """Import harmonypy lazily so the module loads even when the dep is
    absent (e.g. on a slim install)."""
    try:
        import harmonypy
        return harmonypy
    except ImportError as exc:  # pragma: no cover - exercised via runtime check
        raise ImportError(
            "harmonypy is required for Harmony batch correction but is not "
            "installed. It is intentionally removed from the default [bio] "
            "extra because its native build fails on a fresh Windows install. "
            "Install it manually if you need batch correction:\n"
            "  pip install harmonypy"
        ) from exc


def run_harmony(
    adata: AnnData,
    batch_key: str,
    *,
    basis: str = "X_pca",
    corrected_basis: str = "X_pca_harmony",
    max_iter_harmony: int = 10,
    theta: Optional[float] = None,
    random_state: int = 0,
) -> AnnData:
    """Run Harmony batch correction on an embedding.

    Parameters
    ----------
    adata
        AnnData with ``obsm[basis]`` populated (e.g. ``X_pca`` after
        :func:`bioforge.omics.cluster.pca`).
    batch_key
        Name of the ``obs`` column whose values identify batches.
    basis
        Name of the ``obsm`` matrix to correct (default ``X_pca``).
    corrected_basis
        Name of the new ``obsm`` slot to store the corrected embedding.
    max_iter_harmony
        Maximum Harmony iterations to run.
    theta
        Diversity clustering penalty. ``None`` lets harmonypy choose
        the default (``2.0`` per batch).
    random_state
        Reproducibility seed.

    Returns
    -------
    AnnData
        The same object with ``obsm[corrected_basis]`` populated.
    """
    if basis not in adata.obsm:
        raise KeyError(
            f"basis '{basis}' not found in adata.obsm. Call "
            f"bioforge.omics.cluster.pca() first."
        )
    if batch_key not in adata.obs:
        raise KeyError(f"batch key '{batch_key}' not found in adata.obs.")

    harmonypy = _import_harmonypy()

    logger.info(
        "running Harmony: basis=%s batch_key=%s cells=%d comps=%d",
        basis, batch_key, adata.n_obs, adata.obsm[basis].shape[1],
    )
    harmony_out = harmonypy.run_harmony(
        adata.obsm[basis],
        adata.obs,
        batch_key,
        max_iter_harmony=max_iter_harmony,
        theta=[theta] * adata.obs[batch_key].nunique() if theta is not None else None,
        random_state=random_state,
        verbose=False,
    )
    # harmonypy 2.0 already returns Z_corr as (n_cells, n_comps), so it is
    # stored as-is; transposing would misalign cells and components.
    adata.obsm[corrected_basis] = np.asarray(harmony_out.Z_corr)
    logger.info(
        "Harmony complete: corrected embedding stored at obsm['%s'] (shape %s)",
        corrected_basis, adata.obsm[corrected_basis].shape,
    )
    return adata


__all__ = ["run_harmony"]
