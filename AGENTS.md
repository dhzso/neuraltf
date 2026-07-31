# NeuralTF — Planarian Neural TF Discovery Pipeline

## Objective

Rank candidate transcription factors (TFs) that drive neural cell fate in planarians
by integrating 3 single-cell atlas datasets with precomputed enrichment data. Filter
to **neural-fate-specific** candidates using the King G0 progenitor TF atlas as ground
truth, and separate already RNAi-validated genes from novel candidates.

## Scientific Rationale

King et al. 2024 demonstrated that neural fate diversification happens overwhelmingly
**post-mitotically**: 7 broad neoblast clusters expand into 77 G0 progenitor clusters
mapping onto 70+ neuron subtypes. We use King's G0 progenitor compartment assignments
(mmc7) as direct evidence of neural-fate specificity — a TF enriched in a `neural*`
subcluster (log2FC >= 2) is a neural TFs.

### Three atlases provide independent confirmation:

| Atlas | Cells | Genome | Role |
|-------|-------|--------|------|
| **Fincher 2018** | 50,562 | dd_Smed_v4 | scRNA-seq atlas, 44 clusters |
| **Plass 2018** | 21,612 | dd_Smed_v6 | Independent replication atlas |
| **King 2024** | G0/X1 FAC-sorted | dd_Smed_v6 | **Neural ground truth**: TF enrichment across 954 G0 subclusters |

Signal found in all 3 atlases -> high reproducibility. King provides the strongest neural signal but Fincher/Plass confirm the expression pattern is real across independent experiments and closterings.

## Data Sources

| Source | Type | Format | Location |
|--------|------|--------|----------|
| **Fincher atlas** (2018) | scRNA-seq, 26K genes x 33K cells | h5ad (v4 gene IDs) | `datasets/processed/fincher_subsample.h5ad` |
| **Plass atlas** (2018) | scRNA-seq, 17K genes x 37K cells | h5ad (v6 gene IDs) | `datasets/processed/plass_v6.h5ad` |
| **King TF catalog** (2024) | mmc4: 418 TFs with GenBank names | xlsx | `datasets/raw/Supplementary_Data_ King_2024/mmc4.xlsx` |
| **King RNAi** (2024) | mmc5: 162 RNAi-tested FSTFs | xlsx | `datasets/raw/Supplementary_Data_ King_2024/mmc5.xlsx` |
| **King TF pair correlations** (2024) | mmc6: G0-X1 correlation gain | xlsx | `datasets/raw/Supplementary_Data_ King_2024/mmc6.xlsx` |
| **King TF Atlas** (2024) | mmc7: G0, X1, X1 Major Tissue enrichment | xlsx -> tsv | `projects/NeuralTF/data/king_atlas.tsv` |
| **Rosetta Stone** | v4<->v6 gene ID mapping | txt | `datasets/raw/smed_20140614.mapping.rosettastone.2020.txt` |
| **Bridge** | v4<->v6<->gene_name lookup for TF genes | csv | `projects/NeuralTF/data/bridge.csv` |

## Pipeline Steps

1. **Load** Fincher + Plass h5ad files (subsample to 10,000 cells per atlas)
2. **Load references**: TF catalog, RNAi table, correlations, bridge (enrich blank gene_names from mmc4)
3. **QC**: filter, normalize, log1p, HVG (force-include TF genes), PCA, leiden clustering
4. **Score per atlas**: DE via `rank_genes_groups` -> expression (log2FC/5) + specificity (1/n_clusters)
5. **King atlas**: seed all TFs enriched in neural G0 subclusters with log2FC >=2,
   compute expression, specificity, neural_enriched, neural_specificity
6. **RNAi**: match candidates against mmc5 using short dd-ID lookup or gene name
7. **Correlations**: match against mmc6 G0-X1 correlation pairs using ID resolution
8. **Reproducibility**: fraction of atlases supporting each candidate
9. **Function (stub)**: hardcoded known neural TF list
10. **Tier**: RNAi-validated = HIGH-, plus integrated score thresholds
11. **Output**: rank.csv (all candidates), rank_neural.csv (neural-filtered), evidence_cards.md, pipeline_results.json

## Evidence Streams (8 tiers)

| Stage | Weight | Score Computation |
|-------|--------|-------------------|
| Expression | 0.20 | max across all atlases (log2FC/5 in Fincher/Plass, FC/kingMax in King) |
| Specificity | 0.10 | 1 / normalized count percent |
| Reproducibility | 0.15 | n_atlases_supporting / 3 |
| RNAi | 0.15 | 1 if gene_id or short ID in parsed mmc5 targets |
| Correlation | 0.10 | Gain = G0_corr - X1_corr; score = min(1.0, gain * 3) |
| Function | 0.05 | 1 if gene_name in well-known neural list, 0 otherwise |
| **neural_enriched** | 0.15 | 1 if gene hit in any `neural*` G0 subcluster with log2FC >= 2 |
| **neural_specificity** | 0.10 | 1 / number of unique neural subclusters |

## Tier Assignment

- **High**: RNAi valid (rnai > 0) ***OR*** (streams >= 3 AND score >= 0.45)
- **Medium**: streams >= 2 AND score >= 0.25
- **Low**: streams < following

## How to Run

```bash
.venv\Scripts\activate
python -m bioforge.projects.neuraltf.pipeline
```

## Output Files

All written to `projects/NeuralTF/runs/pipeline_run/`:

| File | Content |
|------|---------|
| `rank_neural.csv` | **Neural-enriched candidates ONLY** with proof_status |
| `rank.csv` | All 160+ candidates ranked by integrated score |
| `evidence_cards.md` | Per-candidate evidence cards |
| `ai_summary.md` | AI summary (stub, needs API config) |
| `pipeline_results.json` | Machine-readable top 50 |

## Proof Statuses

- `known_rnai_validated` — in King RNAi table (mmc5) -> existing phenotype
- `prior_fstf_not_tested` — known FSTF from literature, not RNAi' tested possible
- `novel_candidate` — not in mmc5, not a prior FSTF -> priority for validation

## Key Files

| File | Purpose |
|------|---------|
| `src/bioforge/projects/neuraltf/pipeline.py` | Main pipeline class (NeuralTFPipeline) |
| `src/bioforge/evidence/schema.py` | EvidenceRecord, EvidenceSource, ConfidenceTier |
| `src/bioforge/evidence/scoring.py` | EvidenceScorer (multi-criterion weighting) |
| `src/bioforge/evidence/confidence.py` | ConfidencePolicy + assign_tiers (RNAi-first to HIGH) |
| `src/bioforge/evidence/cards.py` | build_evidence_card, classify_proof_status |
| `src/bioforge/evidence/gene_mapping.py` | BridgeTable (load_bridge, v4/v6 conversion) |
| `projects/NeuralTF/data/king_atlas.tsv` | Prebuilt King atlas from mmc7 |