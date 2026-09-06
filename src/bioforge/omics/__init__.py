"""BioForge omics subpackage — generic single-cell workflows.

Wraps scverse (scanpy/anndata/scVelo/CellRank) and harmonypy behind a stable
BioForge API. All operations mutate the AnnData in place (scanpy convention)
and return the same object for chaining.

Modules
-------
- :mod:`bioforge.omics.qc`         — quality control and filtering
- :mod:`bioforge.omics.normalize`  — normalization and feature selection
- :mod:`bioforge.omics.cluster`    — PCA, neighborhood, leiden, UMAP
- :mod:`bioforge.omics.trajectory` — PAGA, RNA velocity, CellRank helpers
- :mod:`bioforge.omics.batch`      — Harmony correction (scanpy-bug-free)
"""
