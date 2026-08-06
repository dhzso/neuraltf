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

# 4. Annotate + prioritize RNAi targets (after 3)
python scripts/query_planmine.py             # PlanMine GO/domains/homology (needs internet)
python scripts/prioritize_neural_tfs.py      # dual-track top-10 -> projects/NeuralTF/results/
```

Pipeline outputs are written to `projects/NeuralTF/runs/pipeline_run/`:

| File | Content |
|------|---------|
| `rank_neural.csv`       | Neural-enriched candidates with proof_status |
| `rank.csv`              | All 160+ candidates |
| `evidence_cards.md`     | Per-candidate evidence summary |
| `pipeline_results.json` | Machine-readable top 50 |

### Prioritization outputs (RNAi shortlist)

Run after the pipeline (step 4 above). Committed to `projects/NeuralTF/results/`:

| File | Content |
|------|---------|
| `top10_neural_tfs_prioritized.csv` | 5 Track A + 5 Track B shortlist |
| `candidate_summary_report.md`      | Per-candidate evidence + wet-lab suggestion |

Supporting data committed to `datasets/processed/`:
`planmine_annotations.parquet` (long-format GO/domain/BLAST rows) and
`planmine_transcripts.fasta` (transcript sequences, one record per candidate
that has one). Both are outputs of `scripts/query_planmine.py`.

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

## Dual-track prioritization (top-10 RNAi shortlist)

After the pipeline ranks ~96 neural-enriched candidates, two scripts turn that
rank list into a concrete wet-lab shortlist:

```bash
python scripts/query_planmine.py          # 1. Fetch annotations from PlanMine
python scripts/prioritize_neural_tfs.py   # 2. Score + split into tracks
```

**Step 1 — PlanMine annotation** (`scripts/query_planmine.py`) queries the
PlanMine InterMine API for every `dd_Smed_v6_*` candidate and stores GO terms,
protein domains (PFAM/InterPro), cross-species BLAST hits and the transcript
sequence in `datasets/processed/planmine_annotations.parquet` +
`planmine_transcripts.fasta`. Internal implementation:
`src/bioforge/projects/neuraltf/planmine.py`. Requires internet access to
`planmine.mpibpc.mpg.de` (the committed parquet/FASTA are cached, so the
prioritization step works offline).

**Step 2 — transparent dual-track scoring** (`scripts/prioritize_neural_tfs.py`):

- merges the mmc4 TF catalog (gene symbol, human ortholog, known-TF flag) and
  the PlanMine annotations onto the 96 candidates;
- maps each gene to its v4 alias via `projects/NeuralTF/data/bridge.csv`,
  flagging `unique` / `ambiguous` mappings (no numeric-ID guessing);
- applies small, additive, documented bonuses: TF domain `+0.05`, neural GO
  `+0.03`, TF GO `+0.02`, human ortholog `+0.02`, RNAi-validated `+0.02`
  (capped at 1.0, each category counted once);
- **Track A** = top 5 RNAi-validated candidates; **Track B** = top-5 novel
  candidates that pass a *tangible TF identity* filter (a DNA-binding-domain
  hit in PlanMine or an mmc4 TF flag — no hypothetical factors without a
  domain);
- appends cross-stage dynamics: Plass X1 neoblast mean vs G0 progenitor
  log2FC (requires `plass_v6.h5ad`; add `--skip-x1` to omit).

Scoring rules live in `src/bioforge/projects/neuraltf/prioritize.py`, which are unit-tested. The Streamlit UI shows the same tables under the
**Prioritization** page (`http://localhost:8501`).

---

## Datasets (scripts only — raw data not bundled)

Raw GEO downloads are not committed (they sum to ~2 GB). Build the
processed h5ad files locally:

**Fincher (GSE111764):**

1. Download from <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111764>
2. Extract `GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz` to
   `datasets/raw/GSE111764_GEO_Fincher_atlas/`
3. `python scripts/convert_fincher.py`

**Plass (GSE103633):**

1. Download the `GSE103633_RAW.tar` supplementary file from <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103633>
2. Place it anywhere under `datasets/raw/` (e.g. `datasets/raw/GSE103633_GEO_Plass_atlas/`); the script auto-locates the tar
3. `python scripts/consolidate_plass.py`

**King 2024 supplementary (Cell Reports):**

1. Download mmc4-mmc7.xlsx from the Cell Reports paper supplementary
2. Place under `datasets/raw/Supplementary_Data_ King_2024/`

**Rosetta Stone gene-ID bridge (v4↔v6):**

1. Download from <https://planosphere.stowers.org/pub/analysis/rosetta/smed_20140614.mapping.rosettastone.2020/smed_20140614.mapping.rosettastone.2020.txt>
2. Place at `datasets/raw/smed_20140614.mapping.rosettastone.2020.txt`
3. `python scripts/build_bridge.py`

The SRA records for the raw reads are Plass: SRP117156 (BioProject
PRJNA403817) and Fincher: SRP135258 (BioProject PRJNA438083).

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
│   │   ├── planmine.py                       PlanMine InterMine REST client + classifiers
│   │   └── prioritize.py                     Dual-track scoring (pure, unit-tested)
│   ├── omics/                                scRNA-seq operations
│   ├── workflow/                             YAML workflow engine
│   ├── ai/                                   AI assistant layer
│   ├── cli/                                  Command-line interface
│   ├── ui/                                   Streamlit UI
│   └── core/                                 Config, datasets, logging, plugins
│
├── datasets/
│   ├── processed/                            Built h5ad files (gitignored) +
│   │                                        planmine_annotations.parquet, .fasta (committed)
│   ├── raw/                                  Raw GEO downloads (gitignored)
│   └── reference/                            Reference tables (gitignored)
│
├── projects/NeuralTF/
│   ├── data/
│   │   ├── bridge.csv                        v4<->v6<->gene_name Rosetta stone
│   │   └── king_atlas.tsv                    Prebuilt G0 enrichment data
│   ├── scripts/visualize_results.py          Generate published figures
│   ├── figures/                              12 generated PNGs
│   ├── results/                              top10_neural_tfs_prioritized.csv + summary report
│   └── runs/                                 Pipeline output directory
│
├── scripts/                                   Utility scripts
│   ├── convert_fincher.py                     Build fincher_subsample.h5ad
│   ├── consolidate_plass.py                   Build plass_v6.h5ad
│   ├── build_bridge.py                        Rebuild bridge.csv
│   ├── build_king_atlas.py                    Rebuild king_atlas.tsv from mmc7
│   ├── audit_king_atlas.py                    Diagnostic: King atlas stats
│   ├── query_planmine.py                      Fetch PlanMine annotation/FASTA for candidates
│   ├── prioritize_neural_tfs.py               Build dual-track top-10 shortlist
│   └── run.py                                Pipeline runner (alternate entry)
│
├── tests/                                     Test suite (195 tests passing)
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
    (pages: Run / Results / Prioritization / Assistant)

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
- Gene IDs are bridged via a mandatory Rosetta Stone CSV
- Subsampling uses `random_state=42` everywhere
- AI operations use a deterministic `StubAssistant` unless an API key is configured
- All 195 unit tests pass on a clean install; 14 more skip when optional deps absent

## License

MIT — see [LICENSE](LICENSE)
