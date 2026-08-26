"""Convert Cui 2023 h5ad from SMED IDs to dd_Smed_v6 for pipeline integration.

Reads the raw Cui h5ad (SMED var_names), maps to v6 via Rosetta Stone,
deduplicates 1-to-many mappings (highest total counts wins), and saves
a v6-named h5ad that the pipeline can load directly.

Output: datasets/processed/cui_v6.h5ad
"""
from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from bioforge.projects.neuraltf.smapping import _load_rosetta

RAW = ROOT / "datasets" / "raw" / "OMIX003867_OMIX_Cui_atlas" / "OMIX003867-01" / "singlecell_h5ad"
OUT = ROOT / "datasets" / "processed" / "cui_v6.h5ad"


def build() -> int:
    h5 = RAW / "adata_scRNA_Annotated.h5ad"
    if not h5.exists():
        print(f"  (missing {h5}, skipping)")
        return 0

    print(f"Loading {h5} ...")
    adata = ad.read_h5ad(h5)
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes (SMED IDs)")

    gene_names = [str(v) for v in adata.var_names]

    # Build SMED -> v6 mapping
    print("  Loading Rosetta Stone mapping ...")
    rosetta = _load_rosetta()
    smed_to_v6: dict[str, list[str]] = {}
    for _, row in rosetta.iterrows():
        s = row["smed_id"]
        v = row["v6_id"]
        smed_to_v6.setdefault(s, []).append(v)
    print(f"  Rosetta Stone: {len(smed_to_v6)} SMED IDs")

    # Map each gene to its v6 ID(s)
    mapped_genes: dict[str, list[int]] = {}  # v6_id -> list of column indices
    unmapped = 0
    for idx, smed_id in enumerate(gene_names):
        v6s = smed_to_v6.get(smed_id, [])
        if not v6s:
            unmapped += 1
            continue
        # Use first v6 ID (1:1 RBH preferred)
        v6 = v6s[0]
        mapped_genes.setdefault(v6, []).append(idx)

    print(f"  Mapped: {len(mapped_genes)} v6 IDs, unmapped: {unmapped}")

    # Deduplicate: for 1-to-many, keep the gene with highest total counts
    X = adata.X
    if hasattr(X, "toarray"):
        import scipy.sparse as sp
        if not sp.issparse(X):
            X = sp.csr_matrix(X)

    keep_cols = []
    keep_names = []
    for v6_id, col_indices in mapped_genes.items():
        if len(col_indices) == 1:
            keep_cols.append(col_indices[0])
            keep_names.append(v6_id)
        else:
            # Pick the column with highest total counts
            totals = np.array([X[:, i].sum() for i in col_indices])
            best = col_indices[int(np.argmax(totals))]
            keep_cols.append(best)
            keep_names.append(v6_id)

    # Subset and rename
    X_sub = X[:, keep_cols]
    new_adata = ad.AnnData(
        X=X_sub,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=keep_names),
    )
    new_adata.var_names.name = None

    OUT.parent.mkdir(parents=True, exist_ok=True)
    new_adata.write_h5ad(OUT)
    print(f"\nSaved {new_adata.shape[1]} genes (v6) x {new_adata.shape[0]} cells to {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(build())
