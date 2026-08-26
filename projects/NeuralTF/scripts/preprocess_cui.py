"""Preprocess Cui 2023 scRNA-seq atlas for NeuralTF pipeline.

Loads adata_scRNA_Annotated.h5ad (55K cells, 61 subtypes), maps SMED IDs
to dd_Smed_v6, and computes per-TF expression/specificity statistics.

Output: data/cui_atlas_summary.csv
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
OUT = ROOT / "projects" / "NeuralTF" / "data" / "cui_atlas_summary.csv"


def build():
    h5 = RAW / "adata_scRNA_Annotated.h5ad"
    print(f"Loading {h5} ...")
    adata = ad.read_h5ad(h5)
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

    gene_names = [str(v) for v in adata.var_names]
    X = adata.X

    # Pre-build SMED -> v6 mapping (single load of Rosetta Stone)
    print("  Loading Rosetta Stone mapping ...")
    rosetta = _load_rosetta()
    smed_to_v6_dict: dict[str, list[str]] = {}
    for _, row in rosetta.iterrows():
        s = row["smed_id"]
        v = row["v6_id"]
        smed_to_v6_dict.setdefault(s, []).append(v)
    print(f"  Rosetta Stone: {len(smed_to_v6_dict)} SMED IDs")

    # Gene mapping: SMED -> v6 (only for genes in the h5ad)
    smed_ids = gene_names
    print(f"  Mapping {len(smed_ids)} SMED IDs to v6 ...")
    v6_map: dict[str, list[str]] = {}
    for sid in smed_ids:
        v6s = smed_to_v6_dict.get(sid, [])
        if v6s:
            v6_map[sid] = v6s
    print(f"  Mapped: {len(v6_map)} SMED IDs -> v6")

    # Only keep genes that map to v6 (use dict for O(1) lookup instead of list.index O(N))
    mapped_smed = list(v6_map.keys())
    name_to_idx = {name: i for i, name in enumerate(gene_names)}
    mapped_idx = [name_to_idx[s] for s in mapped_smed if s in name_to_idx]
    print(f"  Processing {len(mapped_idx)} mapped genes (out of {len(gene_names)})")

    # Subset the sparse matrix to only mapped genes (keep sparse for memory)
    X_sub = X[:, mapped_idx]
    if hasattr(X_sub, "toarray"):
        import scipy.sparse as sp
        if not sp.issparse(X_sub):
            X_sub = sp.csr_matrix(X_sub)
    print(f"  Subset matrix: {X_sub.shape}, nnz: {X_sub.nnz if hasattr(X_sub, 'nnz') else '?'}")

    # Cell type annotations
    cell_types = adata.obs["Annotation"].astype(str).tolist()
    big_types = adata.obs["BigCellType"].astype(str).tolist()
    timepoints = adata.obs["TimePoint"].astype(str).tolist()

    obs_df = pd.DataFrame({
        "cell_type": cell_types,
        "big_type": big_types,
        "timepoint": timepoints,
    })

    unique_big = sorted(set(big_types) - {"nan", ""})
    unique_annot = sorted(set(cell_types) - {"nan", ""})
    neural_annot = [c for c in unique_annot if "neuron" in c.lower() or "neural" in c.lower()]
    unique_tp = sorted(set(timepoints) - {"nan", ""})
    print(f"  Big types: {len(unique_big)}, Annotations: {len(unique_annot)}")
    print(f"  Neural annotations: {neural_annot}")

    # Vectorized computation using sparse matrix operations
    # Precompute boolean masks as sparse row selectors
    import scipy.sparse as sp

    # Build indicator matrices for each big type (n_cells x n_big_types)
    big_type_names = list(unique_big)
    bt_indicator = sp.csr_matrix(
        np.column_stack([(obs_df["big_type"] == bt).values for bt in big_type_names])
    )
    # Neural big type mask (column indices in bt_indicator)
    neural_bt_idx = [i for i, bt in enumerate(big_type_names) if "neur" in bt.lower()]

    # Per-timepoint indicator
    tp_names = list(unique_tp)
    tp_indicator = sp.csr_matrix(
        np.column_stack([(obs_df["timepoint"] == tp).values for tp in tp_names])
    )

    # Per-neural-annotation indicator
    neural_annot_names = list(neural_annot)
    if neural_annot_names:
        na_indicator = sp.csr_matrix(
            np.column_stack([(obs_df["cell_type"] == ca).values for ca in neural_annot_names])
        )
    else:
        na_indicator = None

    # Cell counts per big type and timepoint (for mean computation)
    bt_counts = np.asarray(bt_indicator.sum(axis=0)).flatten().astype(float)
    tp_counts = np.asarray(tp_indicator.sum(axis=0)).flatten().astype(float)
    na_counts = np.asarray(na_indicator.sum(axis=0)).flatten().astype(float) if na_indicator is not None else np.array([])

    # Overall mean per gene: X_sub.T @ ones / n_cells
    n_cells = float(X_sub.shape[0])
    overall_means = np.asarray(X_sub.mean(axis=0)).flatten()

    # Per-big-type means: (X_sub.T @ bt_indicator) / bt_counts  -> (n_genes, n_big_types)
    bt_means = (X_sub.T @ bt_indicator).toarray() / np.maximum(bt_counts, 1)  # shape (n_genes, n_big_types)

    # Per-timepoint means
    tp_means_all = (X_sub.T @ tp_indicator).toarray() / np.maximum(tp_counts, 1)

    # Per-neural-annotation means
    if na_indicator is not None and len(neural_annot_names) > 0:
        na_means_all = (X_sub.T @ na_indicator).toarray() / np.maximum(na_counts, 1)
    else:
        na_means_all = None

    # Positive expression mask (for median computation) — use csr for row slicing
    X_csr = sp.csr_matrix(X_sub)

    records = []
    for local_idx, smed_id in enumerate(mapped_smed):
        v6_ids = v6_map[smed_id]
        overall_mean = float(overall_means[local_idx])

        if overall_mean == 0:
            continue

        # Big-type means for this gene
        bmeans = {bt: float(bt_means[local_idx, i]) for i, bt in enumerate(big_type_names) if bt_counts[i] > 0}
        neural_max = max((bt_means[local_idx, i] for i in neural_bt_idx), default=0)

        # Expression score
        fold_changes = {bt: m / overall_mean for bt, m in bmeans.items() if m > 0}
        max_fc = max(fold_changes.values()) if fold_changes else 0
        expression_score = min(1.0, max_fc / 5.0) if max_fc > 0 else 0.0

        # Specificity
        expr_row = X_csr.getrow(local_idx)
        pos_vals = expr_row.data
        median_expr = float(np.median(pos_vals)) if len(pos_vals) > 0 else 0
        expressed_types = [bt for bt, m in bmeans.items() if m > median_expr]
        specificity_score = 1.0 / len(expressed_types) if expressed_types else 0.0

        # Neural enrichment
        neural_enriched = neural_max > median_expr * 2 if median_expr > 0 else False

        # Neural specificity
        if na_means_all is not None:
            nmeans = {ca: float(na_means_all[local_idx, i]) for i, ca in enumerate(neural_annot_names) if na_counts[i] > 0}
            neural_expressed = [ca for ca, m in nmeans.items() if m > median_expr]
        else:
            neural_expressed = []
        neural_specificity_score = 1.0 / len(neural_expressed) if neural_expressed else 0.0

        # Per-timepoint means
        tp_m = {tp: float(tp_means_all[local_idx, i]) for i, tp in enumerate(tp_names) if tp_counts[i] > 0}

        for v6 in v6_ids:
            records.append({
                "gene_id": v6,
                "smed_id": smed_id,
                "expression_score": round(expression_score, 4),
                "specificity_score": round(specificity_score, 4),
                "neural_enriched": neural_enriched,
                "neural_specificity_score": round(neural_specificity_score, 4),
                "max_fold_change": round(max_fc, 4),
                "n_expressed_types": len(expressed_types),
                "n_neural_expressed": len(neural_expressed),
                "mean_expression": round(overall_mean, 4),
                **{f"mean_{tp}": round(float(tp_m.get(tp, 0)), 4) for tp in tp_names},
            })

    df = pd.DataFrame(records)
    # Keep best v6 per gene (highest expression)
    df = df.sort_values("expression_score", ascending=False).drop_duplicates(subset="gene_id", keep="first")
    df = df.sort_values("gene_id").reset_index(drop=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"\nSaved {len(df)} genes to {OUT}")
    print(f"  Neural enriched: {df['neural_enriched'].sum()}")
    print(f"  Expression score > 0: {(df['expression_score'] > 0).sum()}")


if __name__ == "__main__":
    build()
