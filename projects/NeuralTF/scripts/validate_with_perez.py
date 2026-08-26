"""Validate NeuralTF top-10 candidates against Perez 2025 ANANSE predictions.

Cross-references our prioritized TFs against the ANANSE-predicted
TF-target regulatory network from Perez et al. 2025 (MOESM22).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from bioforge.projects.neuraltf.smapping import v6_to_h1smcg, batch_v6_to_h1smcg

RESULTS = ROOT / "projects" / "NeuralTF" / "results"
MOESM22 = (
    ROOT / "datasets" / "raw" / "Supplementary_Data_ Perez_2025"
    / "41467_2025_65712_MOESM22_ESM.xlsx"
)


def load_ananse_neuron_targets() -> pd.DataFrame:
    """Load ANANSE predicted TF-target interactions for neuron fate."""
    raw = pd.read_excel(MOESM22, sheet_name="target_lists", dtype=str)
    # Filter to neuron fate
    neuron = raw[raw["Fate"].str.lower() == "neuron"].copy()
    neuron.columns = [c.strip() for c in neuron.columns]
    print(f"  ANANSE neuron targets: {len(neuron)} interactions")
    print(f"  Unique TFs: {neuron['TF (gene ID)'].nunique()}")
    print(f"  Unique targets: {neuron['Target gene (gene ID)'].nunique()}")
    return neuron


def validate_top10():
    """Check which top-10 candidates appear in ANANSE neuron network."""
    print("Loading ANANSE neuron targets...")
    ananse = load_ananse_neuron_targets()

    # Build h1SMcG -> v6 reverse map for targets
    ananse_tfs_h1 = set(ananse["TF (gene ID)"].unique())
    ananse_targets_h1 = set(ananse["Target gene (gene ID)"].unique())

    # Load our top-10 results
    for name, path in [
        ("FIXED", RESULTS / "top10_neural_tfs_prioritized.csv"),
        ("CENTERED", RESULTS / "dirichlet_top10_prioritized.csv"),
        ("UNIFORM", RESULTS / "dirichlet_uniform_top10.csv"),
    ]:
        if not path.exists():
            print(f"\n{name}: {path} not found, skipping")
            continue

        df = pd.read_csv(path)
        v6_ids = df["gene_id_v6"].tolist() if "gene_id_v6" in df.columns else df.get("gene_id", pd.Series()).tolist()
        v6_ids = [str(v) for v in v6_ids if str(v) != "nan"]

        # Map v6 -> h1SMcG
        v6_to_h1_map = batch_v6_to_h1smcg(v6_ids)

        print(f"\n{'='*60}")
        print(f"  {name} TOP 10 vs ANANSE neuron network")
        print(f"{'='*60}")

        as_tf = 0
        as_target = 0
        not_found = 0

        for _, row in df.iterrows():
            v6 = str(row.get("gene_id_v6", row.get("gene_id", "")))
            gene = row.get("gene_name", v6)
            h1 = v6_to_h1_map.get(v6)

            if not h1:
                print(f"  {gene:<20}  no h1SMcG mapping")
                not_found += 1
                continue

            is_tf = h1 in ananse_tfs_h1
            is_target = h1 in ananse_targets_h1

            if is_tf:
                # Get TF's targets
                tf_targets = ananse[ananse["TF (gene ID)"] == h1]
                n_targets = len(tf_targets)
                top_target_names = tf_targets["Target (gene symbol)"].head(3).tolist()
                print(f"  {gene:<20}  ANANSE TF ({n_targets} neuron targets: {', '.join(top_target_names)})")
                as_tf += 1
            elif is_target:
                # Find which TFs target this gene
                targeting = ananse[ananse["Target gene (gene ID)"] == h1]
                tf_names = targeting["TF(gene symbol)"].unique()[:3]
                print(f"  {gene:<20}  ANANSE target of: {', '.join(tf_names)}")
                as_target += 1
            else:
                print(f"  {gene:<20}  not in ANANSE neuron network")
                not_found += 1

        total = len(df)
        print(f"\n  Summary: {as_tf}/{total} as TF, {as_target}/{total} as target, {not_found}/{total} not found")


if __name__ == "__main__":
    validate_top10()
