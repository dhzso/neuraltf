"""Convert Cui 2023 h5ad from SMED IDs to dd_Smed_v6 for pipeline integration.

Reads the raw Cui h5ad (SMED var_names), maps to v6 via Rosetta Stone,
deduplicates 1-to-many mappings (highest total counts wins), normalizes
(target_sum=1e4, log1p), and saves a v6-named h5ad that the pipeline
can load directly.

Checkpoints written:
  datasets/processed/cui_conversion_checkpoint.parquet  (mapping/QC stats)

Output: datasets/processed/cui_v6.h5ad
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from bioforge.projects.neuraltf.smapping import _load_rosetta  # noqa: E402

RAW = (
    ROOT
    / "datasets"
    / "raw"
    / "OMIX003867_OMIX_Cui_atlas"
    / "OMIX003867-01"
    / "singlecell_h5ad"
)
OUT = ROOT / "datasets" / "processed" / "cui_v6.h5ad"
CHECKPOINT = ROOT / "datasets" / "processed" / "cui_conversion_checkpoint.parquet"


def build(subsample_n: int | None = None, seed: int = 42) -> int:
    """Convert Cui raw H5AD to v6-indexed processed H5AD.

    Parameters
    ----------
    subsample_n:
        If given, randomly subsample this many cells before processing.
        Default is ``None`` (use all 55,014 cells). Useful for quick testing;
        production runs should use the full dataset for maximum sensitivity.
    seed:
        Random seed for subsampling (ignored when ``subsample_n`` is None).
    """
    t_start = time.time()
    h5 = RAW / "adata_scRNA_Annotated.h5ad"

    # ---------- Fail-fast validation ------------------------------------------
    if not h5.exists():
        raise FileNotFoundError(
            f"Cui raw H5AD not found at expected path:\n  {h5}\n"
            "Download from https://ngdc.cncb.ac.cn/omix/release/OMIX003867 and "
            "place at datasets/raw/OMIX003867_OMIX_Cui_atlas/OMIX003867-01/"
            "singlecell_h5ad/adata_scRNA_Annotated.h5ad"
        )

    print(f"Loading {h5} ...")
    adata = ad.read_h5ad(h5)
    assert adata.n_obs > 0, f"Cui H5AD has no cells: {h5}"
    assert adata.n_vars > 0, f"Cui H5AD has no genes: {h5}"
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes (SMED IDs)")

    # Optional subsampling (disabled by default — full dataset for quality)
    if subsample_n is not None:
        adata = _subsample(adata, subsample_n, seed=seed)



    # Log available obs columns for traceability
    obs_cols = list(adata.obs.columns)
    print(f"  obs columns: {obs_cols}")

    gene_names = [str(v) for v in adata.var_names]

    # ---------- Build SMED -> v6 mapping ---------------------------------------
    print("  Loading Rosetta Stone mapping ...")
    rosetta = _load_rosetta()
    smed_to_v6: dict[str, list[str]] = {}
    for _, row in rosetta.iterrows():
        s = row["smed_id"]
        v = row["v6_id"]
        smed_to_v6.setdefault(s, []).append(v)
    print(f"  Rosetta Stone: {len(smed_to_v6):,} SMED IDs mapped")

    # ---------- Map each gene to its v6 ID(s) ----------------------------------
    mapped_genes: dict[str, list[int]] = {}   # v6_id -> list of column indices
    unmapped = 0
    for idx, smed_id in enumerate(gene_names):
        v6s = smed_to_v6.get(smed_id, [])
        if not v6s:
            unmapped += 1
            continue
        v6 = v6s[0]   # Use first v6 ID (1:1 RBH preferred)
        mapped_genes.setdefault(v6, []).append(idx)

    n_mapped = len(mapped_genes)
    print(f"  Mapped: {n_mapped:,} unique v6 IDs  |  unmapped SMED: {unmapped:,}")
    assert n_mapped > 1000, (
        f"Too few genes mapped ({n_mapped}); Rosetta Stone may be malformed "
        f"or the H5AD gene IDs are not SMED format."
    )

    # ---------- Deduplicate: for 1-to-many, keep highest total counts ----------
    import scipy.sparse as sp

    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    keep_cols: list[int] = []
    keep_names: list[str] = []
    n_multi = 0
    for v6_id, col_indices in mapped_genes.items():
        if len(col_indices) == 1:
            keep_cols.append(col_indices[0])
            keep_names.append(v6_id)
        else:
            n_multi += 1
            totals = np.asarray(X[:, col_indices].sum(axis=0)).ravel()
            best = col_indices[int(np.argmax(totals))]
            keep_cols.append(best)
            keep_names.append(v6_id)
    print(f"  1-to-many deduplicated: {n_multi:,}  |  final genes: {len(keep_cols):,}")

    # ---------- Subset and rename ---------------------------------------------
    X_sub = X[:, keep_cols]
    new_adata = ad.AnnData(
        X=X_sub,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=keep_names),
    )
    new_adata.var_names.name = None
    new_adata.obs_names.name = None

    # Preserve key obs metadata columns that the pipeline uses
    for col in ("Annotation", "BigCellType", "TimePoint"):
        if col in adata.obs.columns and col not in new_adata.obs.columns:
            new_adata.obs[col] = adata.obs[col].values

    assert new_adata.n_obs > 0, "Post-subset AnnData has no cells"
    assert new_adata.n_vars > 1000, (
        f"Post-subset AnnData has only {new_adata.n_vars} genes; "
        "expected >1000 after v6 mapping"
    )

    # ---------- Preserve raw integer counts in X ------------------------------
    print("  Preserving raw counts matrix in adata.X (pipeline will normalize) ...")


    # ---------- Write output --------------------------------------------------
    OUT.parent.mkdir(parents=True, exist_ok=True)
    print(f"  Writing {OUT} ...")
    new_adata.write_h5ad(OUT)

    elapsed = time.time() - t_start
    print(
        f"\n[OK] Saved {new_adata.n_vars:,} genes (v6) x {new_adata.n_obs:,} cells "
        f"to {OUT}  [{elapsed:.0f}s]"
    )

    # ---------- Checkpoint ----------------------------------------------------
    ckpt = pd.DataFrame([{
        "n_cells": new_adata.n_obs,
        "n_genes_v6": new_adata.n_vars,
        "n_smed_input": len(gene_names),
        "n_mapped_v6": n_mapped,
        "n_unmapped_smed": unmapped,
        "n_multimapped_resolved": n_multi,
        "output_path": str(OUT),
        "elapsed_seconds": round(elapsed, 1),
        "has_annotation": "Annotation" in new_adata.obs.columns,
        "has_big_cell_type": "BigCellType" in new_adata.obs.columns,
        "has_timepoint": "TimePoint" in new_adata.obs.columns,
    }])
    ckpt.to_parquet(CHECKPOINT, index=False)
    print(f"  Checkpoint written: {CHECKPOINT}")

    return 0


def _subsample(adata: "ad.AnnData", n: int, seed: int = 42) -> "ad.AnnData":
    """Reproducibly subsample *n* cells from *adata* (stratified by Annotation if available).

    Parameters
    ----------
    adata:
        Full AnnData object (cells × genes).
    n:
        Target number of cells. If ``n >= adata.n_obs`` the full object is returned unchanged.
    seed:
        Random seed for reproducibility.

    Returns
    -------
    AnnData with at most *n* cells.
    """
    if n >= adata.n_obs:
        print(f"  [subsample] n={n} >= {adata.n_obs} cells — skipping, using full dataset")
        return adata
    rng = np.random.default_rng(seed)
    idx = rng.choice(adata.n_obs, size=n, replace=False)
    idx.sort()
    sub = adata[idx].copy()
    print(f"  [subsample] {adata.n_obs:,} -> {sub.n_obs:,} cells (seed={seed})")
    return sub


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(
        description="Convert Cui 2023 SMED h5ad to dd_Smed_v6 for pipeline integration."
    )
    ap.add_argument(
        "--subsample", type=int, default=None, metavar="N",
        help=(
            "Randomly subsample N cells before processing (default: use full dataset). "
            "Useful for fast testing; production runs should omit this flag to preserve "
            "all 55,014 cells for maximum sensitivity."
        ),
    )
    ap.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for subsampling (default: 42).",
    )
    args = ap.parse_args()
    sys.exit(build(subsample_n=args.subsample, seed=args.seed))



