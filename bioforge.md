# bioforge.md — Operations, Architecture & Deep Dive

Extended reference for BioForge. Quick start: see [README.md](README.md).

---

## 1. Complete operational workflow

### 1.1 First-time setup

```bash
git clone https://github.com/dhzso/neuraltf.git
cd neuraltf
python -m venv .venv
.venv\Scripts\activate        # Windows
# source .venv/bin/activate   # Linux/Mac
pip install -e ".[bio,streamlit]"
```

### 1.2 Verify installation (instant, no data needed)

```bash
python -c "from bioforge import *; print('BioForge OK')"

python -m pytest tests/unit/test_evidence.py tests/unit/test_evidence_cards.py -v -q
python -m pytest tests/unit/test_neuraltf_prioritize.py -q   # planmine+prioritize logic
```

### 1.3 Build datasets from raw GEO downloads

### 1.3 Build datasets from raw downloads

The processed h5ad and reference files are not committed. Build them
locally from raw downloads (see `datasets/MANIFEST.md`):

```bash
python scripts/convert_fincher.py                     # Fincher atlas (GSE111764)
python scripts/consolidate_plass.py                   # Plass atlas (GSE103633)
python projects/NeuralTF/scripts/convert_cui.py       # Cui atlas (OMIX003867, full 55K cells)
python projects/NeuralTF/scripts/preprocess_perez.py  # Perez 2025 TF classes (MOESM5)
python scripts/build_master_catalog.py                # Merge King mmc4 + Perez MOESM5
```

The King 2024 supplementary xlsx files (mmc4–mmc7) and Perez 2025 files (MOESM5, MOESM22)
are auto-discovered under their respective directories in `datasets/raw/`.

### 1.4 Run the NeuralTF pipeline (CLI)

```bash
bioforge neuraltf run [--subsample 0] [--out projects/NeuralTF/runs/my_run]
# Or directly:
python scripts/run.py
```

This runs the full pipeline across all 5 atlases and writes results to `projects/NeuralTF/runs/pipeline_run/`:

| File | Content |
|------|---------|
| `rank.csv` | All **289 TF candidates** ranked by 8-stream integrated score |
| `rank_neural.csv` | **102 neural-enriched candidates** with proof_status |
| `evidence_cards.md` | Per-candidate markdown evidence summary (289 cards) |
| `pipeline_results.json` | Machine-readable top 50 with tier metadata |
| `checkpoint_01` – `06` | Parquet audit checkpoints at each major pipeline step |

### 1.5 Evidence Streams & Scoring Model (8 Streams)

Scoring integrates 8 distinct biological streams with default weights:
- **Expression ($w_1 = 0.200$)**: $\min(1.0, \max(\text{log}_2\text{FC})/5)$ across Fincher, Plass, and Cui scRNA-seq atlases.
- **Specificity ($w_2 = 0.100$)**: $1 / n_{\text{clusters}}$ supporting differential expression.
- **Reproducibility ($w_3 = 0.100$)**: $n_{\text{atlases supporting}} / 3$ (Fincher, Plass, Cui).
- **RNAi ($w_4 = 0.100$)**: 1.0 if functional phenotype observed in King mmc5 screen.
- **Correlation ($w_5 = 0.100$)**: $\min(1.0, \Delta r_{\text{G0-X1}} \times 3.0)$ co-expression gain from King mmc6.
- **Neural Enriched ($w_6 = 0.100$)**: 1.0 for G0 neural subcluster log₂FC ≥ 2.0 in King mmc7.
- **Neural Specificity ($w_7 = 0.100$)**: $1 / n_{\text{neural subclusters}}$ present in King atlas.
- **Perez Lineage ($w_8 = 0.100$)**: Perez 2025 lineage class (1.0 neural-class, 0.5 other TF class, 0.0 absent).

Weights renormalize over present streams per candidate.

### 1.5b Dirichlet-robust prioritization (weight sensitivity)

To test whether the ranking is robust to weight assumptions, we perform Monte Carlo Dirichlet sampling (1,000 draws, seed=2024):

```bash
# Centered Dirichlet (k=40)
python projects/NeuralTF/scripts/dirichlet_centered_all249.py   # All 289 candidates
python projects/NeuralTF/scripts/dirichlet_prioritize.py         # 102 neural candidates

# Uniform Dirichlet (α=1)
python projects/NeuralTF/scripts/dirichlet_uniform_all249.py    # All 289 candidates
python projects/NeuralTF/scripts/dirichlet_uniform.py           # 102 neural candidates
```

**Key Finding:** 10/10 top-10 overlap across fixed-weight, centered Dirichlet, and uniform Dirichlet methods confirms that candidate ranking is exceptionally stable.

### 1.5c ANANSE Gene Regulatory Network Scan

Validates candidate TFs against the Perez 2025 computational GRN across 9 cell fate lineages (13,746 interactions):

```bash
python projects/NeuralTF/scripts/ananse_full_scan.py
```

Outputs:
- `projects/NeuralTF/results/ananse_network_full.csv` (289 candidates; 30 TF regulators, 50 target genes)
- `projects/NeuralTF/results/ananse_top_regulators.csv` (top regulators ranked by network out-degree)

### 1.6 Launch the Streamlit UI

```bash
bioforge ui [--port 8501] [--host localhost]
```

Opens http://localhost:8501. Pages: **Run** (dataset status + live pipeline
run), **Results** (rank CSVs + interactive visualization panel),
**Prioritization** (coverage metrics + Track A/B tables + full markdown
report), **Assistant** (interactive BioForge assistant).

### 1.7 Generate Publication Figures (21 Figures)

```bash
python projects/NeuralTF/scripts/generate_publication_figures.py
# Or run all downstream steps in one go:
python scripts/run_downstream.py
```

Outputs 21 publication-quality PNG figures into `projects/NeuralTF/figures/`.


---

## 2. Building data from raw sources (Single-line Command)

**Nothing is committed — everything is generated locally** from the raw
downloads in `datasets/raw/` (only the *sources* are needed: GEO/SRA
downloads, Rosetta Stone table, go.obo, King xlsx; see `datasets/MANIFEST.md`).
The one-command path is:

```bash
python scripts/generate_all.py            # everything, incl. PlanMine (network)
```

which runs the sections below in dependency order (2.1, 2.2, pipeline,
2.4, prioritization, figures), each gated on its inputs. The sections below
document each build step individually.


### 2.1 Bridge CSV (v4 <-> v6 gene IDs)

```bash
python scripts/build_bridge.py \
  --rosetta datasets/raw/smed_20140614.mapping.rosettastone.2020.txt \
  --mmc4 "datasets/raw/Supplementary_Data_ King_2024/1-s2.0-S2211124724001712-mmc4.xlsx" \
  --out projects/NeuralTF/data/bridge.csv
```

### 2.2 King atlas TSV (from mmc7.xlsx)

```bash
python scripts/build_king_atlas.py \
  --mmc7 "datasets/raw/Supplementary_Data_ King_2024/mmc7.xlsx" \
  --out projects/NeuralTF/data/king_atlas.tsv
```

### 2.3 Raw download inventory

| Dataset | GEO Accession | Files needed |
|---------|--------------|--------------|
| Fincher 2018 | GSE111764 | `PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz` |
| Plass 2018 | GSE103633 | `RAW.tar` (per-sample DGE tar archive, named `GSE103633_RAW.tar` on GEO) |
| King 2024 | Cell Reports | 4 supplementary xlsx: mmc4 through mmc7 |

The Plass per-sample DGE files are consolidated in-memory by
`scripts/consolidate_plass.py` — you only need `RAW.tar` from GEO (accession
GSE103633).

### 2.4 PlanMine annotations (online, cached)

No download needed — fetched live from the PlanMine API:

```bash
python scripts/query_planmine.py          # full run (99 candidates, ~1-2 min)
python scripts/query_planmine.py --limit 5   # smoke test (no rate concerns)
```

Re-run any time; it refreshes the two local (gitignored) files:
`datasets/processed/planmine_annotations.parquet` and
`datasets/processed/planmine_transcripts.fasta`.

---

## 3. AI assistant configuration

The AI Assistant page (`Streamlit` tab) and the `bioforge.ai` module depend
on a **compatible API provider** (OpenAI, OpenRouter, Together, Groq, vLLM,
SGLang, LM Studio, Ollama).

### 3.1 Enable with environment variables

```bash
# Minimal -- leave only the provider key
export BIOFORGE_HOME_KEY="sk-proj-..."
# Full configuration (all optional except the first)
export BIOFORGE_AI_BASE_URL="https://api.openai.com/v1"           # default
export BIOFORGE_AI_MODEL="gn-4o-mini"                             # default
export BIOFORGE_AI_TEMPERATURE="0.0"                               # default
export BIOFORGE_AI_MAX_TOKENS="1024"                               # default
```

### 3.2 Providers known to work

| Provider | Base_UR | Model env that'll work |
|----------|---------|----------------|
| OpenAI API | `https://api.openai.com/v1` | `gpt-4o-mini` default, `gpt-4o`, `gpt-4.1` work |
| OpenRouter | `https://openrouter.ai/api/v1` | `openai/gpt-4o-mini`, `google/gemini-2.0-flash` |
| Ollama (local) | `https://localhost:11434/v1` | `llama3.2:3b`, `deepseek-r1:8b` |
| LM Studio (local) | `https://localhost:1234/v1` | Whatever model you load |

### 3.3 Striking the assistant without a key (stub mode)

If no key is set, the assistant falls back to a deterministic `StubAssistant`
that returns canned responses. This is intentional — the pipeline and UI
both work out of the box with no API key. Only the quality of AI-generated
summaries degrades.

### 3.4 Testing connectivity

```bash
python -c "
from bioforge.ai import build_assistant
from bioforge.ai.assistant import ChatMessage
import os; os.environ.setdefault('BIOFORGE_AI_API_KEY', 'test-key')
assistant = build_assistant()
print('Active:', assistant.name)
"
```

---

## 4. Architecture

### 4.1 Layer stack

```
Layer 0-2   Foundational: config, logging, exceptions, plugins
Layer 3     Sci stack: numpy, scipy, pandas
Layer 4     Bioinformatics: scanpy, anndata, leidenalg, igraph
Layer 5     CLI: click-based command interface
Layer 6     AI: OpenAI-compatible provider (stub by default)
Layer 7     Workflow: YAML-engine with step registry + provenance
Layer 8A    Omics: QC, normalize, cluster, trajectory, batch-correction
Layer 8B    Evidence: schema, scoring, cards, gene_mapping, bridge
Layer 8C    Ingest: auto-detect and load h5ad, DGE, TSV, 10X-mtx
Layer 9     NeuralTF: pipeline.py + planmine.py + prioritize.py
Layer 10    UI: Streamlit app (Run / Results / Prioritization / Assistant)
```

### 4.2 Evidence Integration Framework (8B)

```
EvidenceSource (Enum):
  EXPRESSION, SPECIFICITY, REPRODUCIBILITY,
  RNai, CORRELATION,
  NEURAL_ENRICHED, NEURAL_SPECIFICITY

EvidenceRecord                          @dataclass:
  gene_id: str
  gene_name: Optional[str]
  scores: dict[EvidenceSource -> float]
  notes: dict[EvidenceSource -> str]
  proof_status: Optional[str]

EvidenceScorer:
  integrated_score(record) -> float     # weighted renormalized sum

BridgeTable (gene_mapping.py):
  adata.var_names {dd_smed_v6_...} <--> dd_Smed_v6_... {gene symbol}
```

### 4.3 Pipeline call sequence

```
NeuralTFPipeline.run()
  1. load_datasets()                       -> adata_fincher, adata_plass
  2. load_reference_tables()               -> tf_catalog, rnai_table,
                                              correlations, bridge
  3. run_qc()                              -> log1p + HVG + PCA + leiden
  4. score_atlases()                       -> percluster wilcoxon DE, expr + spec
  5. integrate_king_atlas()                -> seed neural-enriched IDs from G0 atlas
  6. integrate_rnai()                      -> match RNAi targets with short dd-ID resolution
  7. integrate_correlations()              -> correlate G0/X1 TF-pairs from mmc6
  8. assign_reproducibility()              -> n_atlases_supporting / 3
  9. write_outputs()                       -> CSV + JSON + cards + terminal panel
```

### 4.4 Prioritization call sequence

```
scripts/query_planmine.py                    (network, one-off; parquet/fasta cached locally)
  1. load_candidates(rank_neural.csv)        -> [gene_id, gene_name]
  2. PlanMineClient.fetch_contig_annotations -> GO / domains / BLAST / sequence
  3. write parquet (long rows) + FASTA

scripts/prioritize_neural_tfs.py
  1. prepare_candidates(rank, mmc4)          -> name/ortholog/TF flag merged
  2. merge_annotations(parquet)              -> PlanMine wins over defaults
  3. map_v4(bridge)                          -> unique / ambiguous / unmapped
  4. compute_composite()                     -> additive bonuses, cap 1.0
  5. assign_tracks()                         -> Track A (RNAi top-5) / B (novel top-5)
  6. X1 vs G0 dynamics (Plass; --skip-x1)    -> cross-stage log2FC + X1 means
  7. write top10 CSV + summary report MD
```

---

## 5. Evidence streams and weights

| Stream | Weight | Interpretation |
|--------|--------|----------------|
| Expression | 0.211 | max log2FC/5 across Fincher and Plass atlases; /8.77 (max across King atlas) |
| Specificity | 0.105 | 1 / n_clusters supporting the TF in this atlas |
| Reproducibility | 0.158 | Fraction (n_atlases_supporting / 3) |
| RNAi | 0.158 | 1 if present in King mmc5 RNAi table |
| Correlation | 0.105 | G0/X1 TF-pair gain capped |
| Neural Enriched | 0.158 | 1 if neural with G0 subcluster log2FC >= 2 |
| Neural Specificity | 0.105 | 1 / n_unique_neural_subclusters (higher = more specific) |

Weights sum to 1.0 (the old `function` stream was removed and its 0.05 was
re-allotted proportionally across the seven streams). The scorer
renormalizes over the streams present per candidate, so the absolute sum is
cosmetic. Note: `neural_enriched` is the cohort-defining neural gate
restated — constant 1.0 within the neural subset; it separates neural from
non-neural candidates in the full 249 list but does not discriminate inside
the 97.

### 5.1 Tier assignment

- **HIGH**: RNAi-validated OR streams > == 3 AND score >= 0.45
- **MEDIUM**: streams >= 2 AND score >= 0.25
- **LOW**: everything else

## 5.2 Proof status

| Status | Meaning | Priority |
|--------|---------|----------|
| `known_rnai_validated` | Already tested in King's RNAi column | Low - known |
| `novel_candidate`     | Neither RNAi nor known FSTF from literature | **HIGH - run experiment** |
| `prior_fstf_not_tested` | Known FSTF, not in RNAxi table | Literature Known|

## 6. Troubleshooting + known issues

### 6.1 `harmonypy` / `scvelo` / `cellrank` not installed

These heavy optional deps are **intentionally** excluded from the default
`[bio]` extra because their native builds require BLAS/CMake and often fail
on a fresh Windows install. If you need batch-correction or trajectory:
```bash
pip install harmonypy  # for bioforge.omics.batch.run_harmony
pip install scvelo     # for bioforge.omics.trajectory.velocity
pip install cellrank  # for bioforge.omics.trajectory.cellrank_terminal_states
```

Without them, the modules import safely, and calling the wrapped functions
raises the correct `ImportError` telling you how to install them.

### 6.2 `MemoryError` during `rank_genes_groups`

`scanpy.tl.rank_genes_groups` uses the full dense matrix which can require
~2 GB with 10K cell x 50K gene matrices. If you hit OOM:
- Reduce `--highly_variable_genes` (pipeline defaults to 5000)
- Use a machine with >= 8 GB RAM
- The default `--subsample 0` keeps the complete atlases; set `--subsample 5000` for lower memory usage / faster development runs

### 6.3 `No such file datasets/processed/fincher_subsample.h5ad`

The h5ad files are not committed. Build them:
```bash
python scripts/convert_fincher.py
python scripts/consolidate_plass.py
```

Make sure raw downloads are in `datasets/raw/` (see section 2).

---

## 9. Directory map

| Path | What |
|------|------|
| `README.md` | Quick start |
| `bioforge.md` | This file — operations + architecture + deep dive |
| `pyproject.toml` | Package config, deps, CLI entry points |
| `src/bioforge/evidence/` | Scoring engine (7 streams), confidence, cards |
| `src/bioforge/projects/neuraltf/pipeline.py` | Main pipeline |
| `src/bioforge/projects/neuraltf/planmine.py` | PlanMine InterMine client + DBD/GO classifiers |
| `src/bioforge/projects/neuraltf/prioritize.py` | Dual-track scoring, mapping, track assignment |
| `src/bioforge/cli/`     | CLI commands (neuraltf, ui, info, run, etc.) |
| `src/bioforge/ui/` | Streamlit app |
| `src/bioforge/ai/` | AI assistant providers |
| `src/bioforge/omics/` | ScRNA-seq operations (QC, normalize, leiden, harmony, etc.) |
| `src/bioforge/workflow/` | Declarative YAML workflow engine |
| `scripts/` | Utility scripts: build_bridge.py, convert_fincher.py, query_planmine.py, prioritize_neural_tfs.py, generate_all.py, etc. |
| `projects/NeuralTF/data/` | bridge.csv, king_atlas.tsv (generated, gitignored) |
| `projects/NeuralTF/results/` | top10_neural_tfs_prioritized.csv + candidate_summary_report.md (generated, gitignored) |
| `datasets/processed/` | h5ad + planmine_annotations.parquet / planmine_transcripts.fasta (generated, gitignored) |
| `projects/NeuralTF/figures/` | 12 visualization figures (generated, gitignored) |
| `projects/NeuralTF/figures/supplementary/` | 4 GO supplementary figures + matrix CSV (generated, gitignored) |
| `projects/NeuralTF/runs/` | Pipeline output dirs (not committed) |
| `projects/NeuralTF/scripts/visualize_results.py` | Figure generator |
| `tests/` | Test suite (195 passing, 14 skipped when optional deps absent) |
| `datasets/` | Raw + processed data (not committed) |
| `.streamlit/config.toml` | No-email config for Streamlit |

---

## 7. Deleting the entire installation

### Remove the docker setup (if you used one)

```bash
docker compose down -v        # removes networks, volumes
docker rmi bioforge-dev       # remove the built image
```

### Remove the package and venv

```bash
# Deactivate if active
python -m venv cleanup --copies
deactivate

pip uninstall bioforge -y
rm -rf .venv      # or deactivate and Remove-Item -Recurse -Force .venv
rm -rf dist build *.egg-info   # built artifacts
```

### Remove the repo directory entirely

```bash
# From outside the repo directory:
rm -rf neuraltf            # Linux/Mac
Remove-Item -Recurse -Force .\neuraltf    # Windows PowerShell
```

That's the complete removal. No hidden daemons or background services —
streamlit stops on Ctrl+C. No data leaks outside the repo root.

---
