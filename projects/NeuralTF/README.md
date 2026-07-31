# NeuralTF — Planarian Neural-Fate Transcription Factor Discovery

The First Consumer for BioForge. Integrates 3 independent planarian scRNA-seq
atlases to identify high-confidence, neural-fate-specific transcription factor
candidates for RNAi validation.

## Layout

```
projects/NeuralTF/
  README.md                 # this file
  data/
    bridge.csv              # v4<->v6<->gene_name bridge table
    king_atlas.tsv          # Prebuilt King mmc7 G0 enrichment data
  scripts/
    visualize_results.py    # Generate 12 publication figures
  figures/                  # Generated figures
  runs/
    pipeline_run/           # Default output directory
      rank.csv              # All TF candidates
      rank_neural.csv       # Neural-enriched candidates only
      evidence_cards.md     # Per-candidate evidence summary
      pipeline_results.json # Machine-readable top 50
```

## Running the pipeline

Requires processed h5ad files and the bridge CSV + king_atlas.tsv (already in this repo).

```bash
python -m bioforge.projects.neuraltf.pipeline
```

Outputs written to `projects/NeuralTF/runs/pipeline_run/`.

### Generate figures from results

```bash
python projects/NeuralTF/scripts/visualize_results.py
```

Outputs 12 PNGs to `projects/NeuralTF/figures/`.

## How it works

The pipeline integrates **8 evidence streams** per TF candidate:

| Stream | Weight | Source |
|--------|--------|--------|
| Expression | 0.20 | Max log2FC/5 across Fincher, Plass, King atlases |
| Specificity | 0.10 | 1 / n_clusters supporting the candidate |
| Reproducibility | 0.15 | Fraction of 3 atlases confirming the TF |
| RNAi | 0.15 | 1 if in King mmc5 RNAi phenotype tablember |
| Correlation | 0.10 | Gmapped G0-X1 gain from King mmc6 |
| Function | 0.05 | Hardcoded known neural TF list |
| Neural Enriched | 0.15 | Binary 1.0 if TF hits a neural G0 subcluster with log2FC >= 2 |
| Neural Specificity | 0.10 | 1 / number of unique neutral subclusters |

Tier assignment:
- **HIGH**: RNAi-validated OR (streams >= 3 AND integrated >= 0.45)
- **MEDIUM**: streams >= 2 AND integrated >= 0.25
- **LOW**: otherwise

Proof status:
- `known_rnai_validated` — in King mmc5 RNAi table (existing phenotype)
- `prior_fstf_not_tested` — known FSTF from literature, not in mmc5
- `novel_candidate` — not in mmc5, not a prior FSTF — priority for validation

## Reproducibility

All source atlases are independent experiments:
- Fincher 2018 (dd_Smed_v4): 50K cell Drop-seq atlas
- Plass 2018 (dd_Smed_v6): 21K cell Drop-seq atlas, independent lab
- King 2024 (dd_Smed_v6): FAC-sorted G0 porter TF atlas with subcluster annotations

All anthology gene IDs bridged via Rosetta Stone mapping; no guessing by numeric prefix.