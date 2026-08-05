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
```

### 1.3 Build datasets from raw GEO downloads

The processed h5ad files are not committed (they sum to ~80 MB). Build them
locally from raw GEO downloads:

```bash
python scripts/convert_fincher.py      # Fincher atlas (GSE111764)
python scripts/consolidate_plass.py    # Plass atlas (GSE109226)
```

The King 2024 supplementary xlsx files (mmc4-mmc7) must be placed under
`datasets/raw/Supplementary_Data_ King_2024/`. The pipeline auto-discovers
them — the exact Elsevier filename (`1-s2.0-S2211124724001712-mmcN.xlsx`) is
tried first, then any `*mmcN.xlsx` in the directory.

### 1.4 Run the NeuralTF pipeline (CLI)

```bash
bioforge neuraltf run [--subsample 10000] [--out projects/NeuralTF/runs/my_run]
```

This runs the full pipeline and writes results:

| File | Content |
|------|---------|
| `rank.csv` | All ~160 TF candidates ranked by integrated score |
| `rank_neural.csv` | Neural-enriched candidate subset with proof_status |
| `evidence_cards.md` | Per-candidate markdown evidence summary |
| `pipeline_results.json` | Machine-readable top 50 with tier metadata |

### 1.5 Launch the Streamlit UI

```bash
bioforge ui [--port 8501] [--host localhost]
```

Opens http://localhost:8501. The Run page shows dataset status and runs
the pipeline with live progress; the Results page tabs the CSVs and has a
full Visualization panel with 12 matplotlib figures; the Assistant page
provides an AI chatbot (see section on AI below).

### 1.6 Generate publication figures (as PNG)

```bash
python projects/NeuralTF/scripts/visualize_results.py
```

Outputs 12 PNGs to `projects/NeuralTF/figures/`.

---

## 2. Building data from raw sources

The two h5ad files and the King atlas TSV are **prebuilt and committed**.
Rebuild them only if the upstream raw data changes.

### 2.1 Bridge CSV (v4 <-> v6 gene IDs)

```bash
python scripts/build_bridge.py \
  --rosetta datasets/raw/smed_20140614.mapping.rosettastone.2020.txt \
  --mmc4 "datasets/raw/Supplementary_Data_ King_2024/mmc4.xlsx" \
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
| Fincher 2018 | GSE111278 | `PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz` |
| Plass 2018 | GSE109226 | `RAW.tar` extracted (per-cell DGE tar archive) |
| King 2024 | Cell Reports | 4 supplementary xlsx: mmc4 through mmc7 |

The Plass per-cell files are consolidated in-memory by
`scripts/consolidate_plass.py` — you only need `RAW.tar` from GEO.

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
Layer 9     NeuralTF: pipeline.py
Layer 10    UI: Streamlit app (Run / Results / Assistant)
```

### 4.2 Evidence Integration Framework (8B)

```
EvidenceSource (Enum):
  EXPRESSION, SPECIFICITY, REPRODUCIBILITY,
  RNai, CORRELATION, FUNCTION,
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

---

## 5. Evidence streams and weights

| Stream | Weight | Interpretation |
|--------|--------|----------------|
| Expression | 0.20 | max log2FC/5 across all 3 atlases |
| Specificity | 0.10 | 1 / n_clusters | supporting the TF in this atlas |
| Reproducibility | 0.15 | Fraction (n_atlases_supporting / 3) |
| RNAi | 0.15 | 1 if present in King mmc5 RNAi table |
| Correlation | 0.10 | G0/X1 TF-pair gain capped |
| Function | 0.05 | 1 if known neural TF (literature) |
| Neural Enriched | 0.15 | 1 if neural with G0 subcluster log2FC >= 2 |
| Neural Specificity | 0.10 | 1 / n_unique_neural_subclusters (higher = more specific) |

### 5.1 Tier assignment

- **HIGH**: RNAi-validated OR streams > == 3 AND score >= 0.45
- **MEDIUM**: streams >= 2 AND score >= 0.25
- **LOW**: everything else

## 5.2 Proof status

| Status | Meaning | Priority |
|--------|---------|----------|
| `known_rnai_validated` | Already tested in King's RNAi column | Low - known |
| `novel_candidate`     | Neither RNAi nor known FSTF from literature | **HIGH - run experiment** |
| `prior_fstf_not_tested` | Known FSTF, not in RNAxi table | Verified - re-test valid |

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
- Set `--subsample 5000` for lower memory usage

### 6.3 `No such file datasets/processed/fincher_subsample.h5ad`

The h5ad files are not committed. Build them:
```bash
python scripts/convert_finch.py
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
| `src/bioforge/evidence/` | Scoring engine (8 streams), confidence, cards |
| `src/bioforge/projects/neuraltf/pipeline.py` | Main pipeline |
| `src/bioforge/cli/`     | CLI commands (neuraltf, ui, info, run, etc.) |
| `src/bioforge/ui/` | Streamlit app |
| `src/bioforge/ai/` | AI assistant providers |
| `src/bioforge/omics/` | ScRNA-seq operations (QC, normalize, leiden, harmony, etc.) |
| `src/bioforge/workflow/` | Declarative YAML workflow engine |
| `scripts/` | Utility scripts: build_bridge.py, convert_fincher.py, etc. |
| `projects/NeuralTF/data/` | bridge.csv, king_atlas.tsv (committed) |
| `projects/NeuralTF/figures/` | 12viz figures (committed) |
| `projects/NeuralTF/runs/` | Pipeline output dirs (not committed) |
| `projects/NeuralTF/scripts/visualize_results.py` | Figure generator |
| `tests/` | Test suite (170 passing, 14 skipped when optional deps absent) |
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

## 8. Reproducibility guarantees

- Explicit bridge table (v4 <-> v6 <-> gene_name) — no numeric-prefix gene ID guessing
- Deterministic random seed per step (seed=42 for subsampling, seed=0 for PCA, etc.)
- Same code + same inputs = same output IDs and same score values
- All 170 unit tests pass on clean install; 14 more tests get skipped when optional deps absent
- Locked wheel file via pip (auto-resolved on install)

---

## 9. What to do after raw data changes

If upstream GEO releases update or King/Cell Reports publishes a corrected xlsx:

1. `Down date / rebuild processed h5ad`:
   ```bash
   python scripts/convert_fincher.py
   python scripts/consolidate_plass.py
   ```
2. Rebuild reference tables (if King 2024 supplementary files change):
   ```bash
   python scripts/build_bridge.py [--rosetta ...] [--out ...]
   python scripts/build_king_atlas.py [--mmc7 ...] [--out ...]
   ```
3. Re-run the pipeline:
   ```bash
   bioforge neuraltf run                # default output dir
   bioforge neuraltf run --out .../run_after_update
   ```
4. Update visualizations:
   ```bash
   python projects/NeuralTF/scripts/visualize_results.py --run ... --out ...
   ```
5. Commit the changed **derived data** (bridge.csv, king_atlas.tsv): no.
   Replaced raw data goes in `datasets/raw/` and processed in `datasets/processed/`;
   these files are git-ignored by .gitignore.

The pipeline is portable: the only thing a new user has to do is download
the raw GEO files into the `datasets/raw/` layout described in the README.