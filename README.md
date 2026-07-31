# BioForge

A reproducible bioinformatics platform for cross-atlas transcription factor
discovery. Built for planarian neural-fate research (NeuralTF thesis project).

## Quick start

```bash
# 1. Set up Python environment
python -m venv .venv
.venv\Scripts\activate        # Windows
pip install -e ".[bio,dev]"

# 2. Run the NeuralTF pipeline
python -m bioforge.projects.neuraltf.pipeline

# 3. See results
#    projects/NeuralTF/runs/pipeline_run/
#      rank.csv              — 160 TF candidates ranked by integrated score
#      rank_neural.csv       — 89 neural-specific candidates with proof status
#      evidence_cards.md     — per-candidate evidence summary
#      pipeline_results.json — top 50 candidates (machine readable)
```

## What it does

Integrates **3 independent planarian single-cell RNA-seq atlases** to identify high-confidence neural-fate-specific transcription factors (TFs):

| Atlas | Year | Cells | Role |
|-------|------|-------|------|
| Fincher | 2018 | 50,562 | Whole-animal cell-type atlas (dd_Smed_v4) |
| Plass | 2018 | 21,612 | Independent replication atlas (dd_Smed_v6) |
| King | 2024 | G0 progenitors | Neural ground truth: enrichment across 954 G0 subclusters |

Each candidate TF receives an integrated score from **8 evidence streams** — expression strength (log2FC), specificity, cross-atlas reproducibility, RNAi phenotype confirmation, neural TF-pair correlation gain, literature function, neural enrichment (King G0 log2FC >= 2.0), and neural subtype specificity.

TFs are marked with a **proof status**:
- `known_rnai_validated` — already tested by King et al. via RNAi
- `novel_candidate` — not yet tested — priority for new experiments
- `prior_fstf_not_tested` — known FSTFs from literature without RNAi data

## Project layout

```
Bioinformatics/
├── pyproject.toml                         # Package config + dependencies
├── README.md                              # This file
├── bioforge.md                            # Extended architecture + operations
│
├── src/bioforge/                          # Python package
│   ├── evidence/                          # 8-stream scoring engine
│   │   ├── schema.py                      # EvidenceRecord, EvidenceSource
│   │   ├── scoring.py                     # Weighted score integration
│   │   ├── confidence.py                  # Tier assignment (HIGH/MEDIUM/LOW)
│   │   ├── cards.py                       # Per-candidate evidence cards
│   │   ├── gene_mapping.py                # v4↔v6 bridge table
│   │   └── readers/                       # Atlas dataset readers
│   ├── projects/neuraltf/pipeline.py      # Main pipeline
│   ├── omics/                             # scRNA-seq operations
│   ├── workflow/                          # YAML workflow engine
│   ├── ai/                                # AI assistant layer
│   ├── cli/                               # Command-line interface
│   └── ui/                                # Streamlit UI
│
├── datasets/
│   ├── processed/                         # Preprocessed h5ad files
│   │   ├── fincher_subsample.h5ad         # Fincher atlas (v4, 10K cells)
│   │   └── plass_v6.h5ad                  # Plass atlas (v6, 37K cells)
│   └── raw/(not tracked)                  # Raw GEO downloads + King xlsx
│
├── projects/NeuralTF/
│   ├── data/
│   │   ├── bridge.csv                     # v4 <-> v6 <-> gene_name
│   │   └── king_atlas.tsv                 # Prebuilt G0 enrichment data
│   ├── scripts/visualize_results.py       # Generate published figures
│   ├── figures/                           # 12 generated PNGs
│   └── runs/pipeline_run/                 # Pipeline output directory
│
├── scripts/                               # Utility scripts
│   ├── build_bridge.py                     # Build bridge CSV
│   ├── build_king_atlas.py                 # Build king_atlas.tsv from mmc7
│   └── audit_king_atlas.py                # Diagnostic: King atlas stats
│
├── tests/                                  # Test suite
│   ├── unit/                               # 18 unit test files
│   └── integration/                        # Integration tests
│
├── docker/                                 # Container setup
├── docs/                                   # Architecture decisions
└── .gitignore
```

## Installation

**Minimum Python 3.11**. Install in a virtual environment:

```bash
python -m venv .venv
# Windows
../venv\Scripts\activate
# Linux/Mac
source . venv/bin/activate

pip install -e ".[bio,dev]"
```

### Running inside Docker (optional)

```bash
docker compose up -d
docker exec bioforge-dev bioforge info
```

## Running the pipeline

### Full real-data pipeline

```bash
python -m bioforge.projects.neuraltf.pipeline
```

Outputs written to `projects/NeuralTF/runs/pipeline_run/`:

| File | Content |
|------|---------|
| `rank_neural.csv` | Neural-enriched candidates only with proof_status |
| `rank.csv` | All 160+ candidates |
| `evidence_cards.md` | Per-candidate evidence summary |
| `ai_summary.md` | AI summary (stub if no API key) |
| `pipeline_results.json` | Machine-readable top 50 |

### Demo pipeline (no raw data needed)

```bash
pipeline run projects/Neural_TF/workflows/demo_pipeline.yaml --out projects/NeuralTF/runs/demo
```

### Regenerate derived data (only if raw files change)

```bash
python scripts/build_bridge.py     # Rebuild via CSV
python scripts/build_king_atlas.py  # Rebuild king_atlas.tsv from mmc7
```

### Generate figures

```bash
python projects/NeuralTF/scripts/visualize_results.py
```

## Running tests

```bash
python -m pytest tests/unit -v -k "evidence or projects"
```

## Evidence streams

| Stream | Weight | How it's computed |
|--------|--------|-------------------|
| Expression | 0.20 | max log2FC/5 across all 3 atlases |
| Specificity | 0.10 | 1 / n_clusters supporting the TF |
| Reproducibility | 0.15 | atlases_supporting / 3 |
| RNAi | 0.15 | 1 if gene is in King mmc5 RNAi table |
| Correlation | 0.10 | G0-X1 pair correlation gain × 3 |
| Function | 0.05 | 1 if known neural TF |
| Neural Enriched | 0.15 | 1 if enrichment for any neural G0 subcluster with log2FC ≥ 2.0 |
| Neural Specificity | 0.10 | 1 / n_neural_subclusters present in |

### Tier assignment

- **HIGH**: RNAi-validated OR (streams ≥ 3 AND score ≥ 0.45)
- **MEDIUM**: streams ≥ 2 AND score ≥ 0.25
- **LOW**: all others

## Reproducibility

- Atlas are independent experiments from the peer-reviewed papers
- Gene IDs are bridged via mandatory Rosetta Stone mapping (no numeric-prefix guessing)
- AI operations use a deterministic stub unless a valid API key is provided
- Full dependency list available in `docker/base/requirements.lock.txt`

## License

MIT — see [LICENSE](LICENSE)