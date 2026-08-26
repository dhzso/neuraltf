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

    # Only keep genes that map to v6
    mapped_smed = list(v6_map.keys())
    mapped_idx = [gene_names.index(s) for s in mapped_smed]
    print(f"  Processing {len(mapped_idx)} mapped genes (out of {len(gene_names)})")

    # Subset the sparse matrix to only mapped genes (much faster)
    X_sub = X[:, mapped_idx]
    if hasattr(X_sub, "toarray"):
        X_sub = X_sub.toarray()
    X_sub = np.asarray(X_sub)
    print(f"  Subset matrix: {X_sub.shape}")

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

    # Build masks for big types and timepoints
    big_masks = {bt: (obs_df["big_type"] == bt).values for bt in unique_big}
    neural_big_masks = {bt: big_masks[bt] for bt in unique_big if "neur" in bt.lower()}
    neural_annot_masks = {ca: (obs_df["cell_type"] == ca).values for ca in neural_annot}
    tp_masks = {tp: (obs_df["timepoint"] == tp).values for tp in unique_tp}

    records = []
    for local_idx, smed_id in enumerate(mapped_smed):
        v6_ids = v6_map[smed_id]
        expr = X_sub[:, local_idx]

        if expr.sum() == 0:
            continue

        overall_mean = expr.mean()

        # Per-big-type mean expression
        big_means = {bt: expr[big_masks[bt]].mean() for bt in unique_big if big_masks[bt].sum() > 0}

        # Expression score
        if overall_mean > 0:
            fold_changes = {bt: m / overall_mean for bt, m in big_means.items() if m > 0}
            max_fc = max(fold_changes.values()) if fold_changes else 0
        else:
            max_fc = 0
        expression_score = min(1.0, max_fc / 5.0) if max_fc > 0 else 0.0

        # Specificity
        pos_expr = expr[expr > 0]
        median_expr = float(np.median(pos_expr)) if len(pos_expr) > 0 else 0
        expressed_types = [bt for bt, m in big_means.items() if m > median_expr]
        specificity_score = 1.0 / len(expressed_types) if expressed_types else 0.0

        # Neural enrichment
        neural_max = max((big_means.get(bt, 0) for bt in neural_big_masks), default=0)
        neural_enriched = neural_max > median_expr * 2 if median_expr > 0 else False

        # Neural specificity
        neural_means = {ca: expr[neural_annot_masks[ca]].mean() for ca in neural_annot if neural_annot_masks[ca].sum() > 0}
        neural_expressed = [ca for ca, m in neural_means.items() if m > median_expr]
        neural_specificity_score = 1.0 / len(neural_expressed) if neural_expressed else 0.0

        # Per-timepoint means
        tp_means = {tp: float(expr[tp_masks[tp]].mean()) for tp in unique_tp if tp_masks[tp].sum() > 0}

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
                "mean_expression": round(float(overall_mean), 4),
                **{f"mean_{tp}": round(float(tp_means.get(tp, 0)), 4) for tp in unique_tp},
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
