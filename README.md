# BioForge · NeuralTF

A reproducible pipeline for **planarian neural-fate-specific transcription factor**
discovery. Integrates three peer-reviewed single-cell RNA-seq atlases
(Fincher 2018, Plass 2018, King 2024) plus the King 2024 RNAi phenotype
screen and neural TF-pair correlation data, scoring each candidate TF on
8 evidence streams and flagging priority targets for RNAi validation.

---

## Quick start

```bash
# 1. Set up Python 3.11+ environment
python -m venv .venv
.venv\Scripts\activate                       # Windows
source .venv/bin/activate                    # Linux/Mac

pip install -e ".[bio,streamlit]"

# 2. Build datasets from raw GEO downloads (Scripts only - see "Datasets" below)
python scripts/convert_fincher.py
python scripts/consolidate_plass.py

# 3. Run the pipeline (CLI or Streamlit UI)
bioforge neuraltf run                        # CLI - fastest, no UI
bioforge ui                                  # Streamlit UI - http://localhost:8501
```

Outputs are written to `projects/NeuralTF/runs/pipeline_run/`:

| File | Content |
|------|---------|
| `rank_neural.csv`       | Neural-enriched candidates with proof_status |
| `rank.csv`              | All 160+ candidates |
| `evidence_cards.md`     | Per-candidate evidence summary |
| `pipeline_results.json` | Machine-readable top 50 |

---

## What it does

The pipeline seeds TF candidates from a King 2024 G0 atlas TF catalog,
scores them on 8 evidence streams, then filters for neural-fate specificity.

### Atlases

| Atlas | Year | Cells | Role |
|-------|------|-------|------|
| Fincher | 2018 | 50,562 | Whole-animal cell-type atlas (dd_Smed_v4) |
| Plass   | 2018 | 21,612 | Independent replication atlas (dd_Smed_v6) |
| King    | 2024 | G0 progenitors | Neural ground truth: enrichment across 955 G0 subclusters |

### Evidence streams (8)

| Stream | Weight | How it's computed |
|--------|--------|-------------------|
| Expression         | 0.20 | max log2FC/5 across all 3 atlases |
| Specificity        | 0.10 | 1 / n_clusters supporting the TF |
| Reproducibility    | 0.15 | atlases_supporting / 3 |
| RNAi               | 0.15 | 1 if gene is in King mmc5 RNAi table |
| Correlation        | 0.10 | G0-X1 pair correlation gain x 3 |
| Function           | 0.05 | 1 if known neural TF |
| Neural Enriched    | 0.15 | 1 if enrichment for any neural G0 subcluster with log2FC >= 2.0 |
| Neural Specificity | 0.10 | 1 / n_neural_subclusters present in |

### Tier assignment

- **HIGH**: RNAi-validated OR (streams >= 3 AND score >= 0.45)
- **MEDIUM**: streams >= 2 AND score >= 0.25
- **LOW**: all others

### Proof status

- `known_rnai_validated` — already tested by King et al. via RNAi
- `novel_candidate` — not yet tested — priority for new experiments
- `prior_fstf_not_tested` — known FSTFs from literature without RNAi data

---

## Datasets (scripts only — raw data not bundled)

Raw GEO downloads are not committed (they sum to ~2 GB). Build the
processed h5ad files locally:

**Fincher (GSE111764):**

1. Download from <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111764>
2. Extract `GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz` to
   `datasets/raw/GSE111764_GEO_Fincher_atlas/`
3. `python scripts/convert_fincher.py`

**Plass (GSE109226):**

1. Download from <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103633>
2. Extract `RAW.tar` to `datasets/raw/Plass_2018/`
3. `python scripts/consolidate_plass.py`

**King 2024 supplementary (Cell Reports):**

1. Download mmc4-mmc7.xlsx from the Cell Reports paper supplementary
2. Place under `datasets/raw/Supplementary_Data_ King_2024/`

---

## Project layout

```
Bioinformatics/
├── pyproject.toml                            Package config + dependencies
├── README.md                                 This file
├── bioforge.md                               Extended operations + architecture
├── .streamlit/config.toml                    Streamlit config (no email prompt)
│
├── src/bioforge/                             Python package
│   ├── evidence/                             8-stream scoring engine
│   │   ├── schema.py                         EvidenceRecord, EvidenceSource
│   │   ├── scoring.py                        Weighted score integration
│   │   ├── confidence.py                     Tier assignment (HIGH/MEDIUM/LOW)
│   │   ├── cards.py                           Per-candidate evidence cards
│   │   ├── gene_mapping.py                   v4<->v6 bridge table
│   │   └── readers/                          Atlas dataset readers
│   ├── projects/neuraltf/pipeline.py         Main pipeline (NeuralTFPipeline)
│   ├── omics/                                scRNA-seq operations
│   ├── workflow/                             YAML workflow engine
│   ├── ai/                                   AI assistant layer
│   ├── cli/                                  Command-line interface
│   ├── ui/                                   Streamlit UI
│   └── core/                                 Config, datasets, logging, plugins
│
├── datasets/
│   ├── processed/                            Built h5ad files (gitignored)
│   ├── raw/                                  Raw GEO downloads (gitignored)
│   └── reference/                            Reference tables (gitignored)
│
├── projects/NeuralTF/
│   ├── data/
│   │   ├── bridge.csv                        v4<->v6<->gene_name Rosetta stone
│   │   └── king_atlas.tsv                    Prebuilt G0 enrichment data
│   ├── scripts/visualize_results.py          Generate published figures
│   ├── figures/                              12 generated PNGs
│   └── runs/                                 Pipeline output directory
│
├── scripts/                                   Utility scripts
│   ├── convert_fincher.py                     Build fincher_subsample.h5ad
│   ├── consolidate_plass.py                   Build plass_v6.h5ad
│   ├── build_bridge.py                        Rebuild bridge.csv
│   ├── build_king_atlas.py                    Rebuild king_atlas.tsv from mmc7
│   ├── audit_king_atlas.py                    Diagnostic: King atlas stats
│   └── run.py                                Pipeline runner (alternate entry)
│
├── tests/                                     Test suite (170 tests passing)
└── docs/                                      Architecture decisions
```

---

## CLI reference

```bash
bioforge --version
bioforge --help

# Main commands
bioforge neuraltf run [--subsample N] [--out DIR]   Run the NeuralTF pipeline
bioforge ui [--port 8501] [--host localhost]       Launch the Streamlit UI

# Supporting commands (kept for advanced use)
bioforge info                                      Show build + config
bioforge datasets list [--category raw|processed]   List datasets
bioforge projects list                              List research projects
bioforge plugins list                               List plugins
bioforge run WORKFLOW.yaml                          Execute a YAML workflow
```

---

## Optional dependencies

The default `[bio]` extra installs the slim stack the pipeline actually
needs: scanpy, anndata, biopython, gseapy, igraph, leidenalg, openpyxl.

Heavy scverse extras (`harmonypy`, `scvelo`, `cellrank`) are *intentionally*
excluded because their native BLAS/CMake builds fail on a fresh Windows
install. The `bioforge.omics.batch` and `bioforge.omics.trajectory` modules
import cleanly without them; calling the wrapper functions raises a clear
`ImportError` if you forgot to install the underlying package.

Install them only if you actually use the corresponding wrapper:

```bash
pip install harmonypy    # bioforge.omics.batch.run_harmony
pip install scvelo       # bioforge.omics.trajectory.velocity
pip install cellrank     # bioforge.omics.trajectory.cellrank_terminal_states
```

---

## Reproducibility

- Atlases are independent experiments from peer-reviewed papers
- Gene IDs are bridged via a mandatory Rosetta Stone CSV (no numeric-prefix guessing)
- Subsampling uses `random_state=42` everywhere
- AI operations use a deterministic `StubAssistant` unless an API key is configured
- All 170 unit tests pass on a clean install; 14 more skip when optional deps absent

## License

MIT — see [LICENSE](LICENSE)
