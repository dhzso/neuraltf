# bioforge.md — Operations, Architecture & Deep Dive

Extended reference for BioForge. Quick start: see [README.md](README.md).

---

## 1. Complete operational workflow

### 1.1 First-time setup

```bash
git clone <repo-url>
cd Bioinformatics
python -m venv .venv
.venv\Scripts\activate        # Windows
# source .venv/bin/activate   # Linux/Mac
pip install -e ".[bio,dev]"
```

### 1.2 Verify installation (instant, no data needed)

```bash
python -c "from bioforge import *; print('BioForge OK')"

python -m pytest tests/unit/test_evidence.py tests/unit/test_evidence_cards.py -v -q
```

### 1.3 Run the NeuralTF pipeline

```bash
python -m bioforge.projects.neuraltf.pipeline
```

This runs in ~5 minutes on a standard laptop and writes results to `projects/NeuralTF/runs/pipeline_run/`:

| File | Content |
|------|---------|
| `rank.csv` | 160 TF candidates ranked by integrated score |
| `rank_neural.csv` | 89 neural-enriched candidates only, with proof_status |
| `evidence_cards.md` | Per-candidate markdown evidence summary |
| `ai_summary.md` | AI-generated summary (stub if no API key) |
| `pipeline_results.json` | Machine-readable top 50 |

### 1.4 View results

```bash
python -c "
import pandas as pd
df = pd.read_csv('projects/NeuralTF/runs/pipeline_run/rank_neural.csv')
print(df[['gene_name','integrated_score','rnai','proof_status']].to_string())
"
```

### 1.5 Generate publication figures

```bash
python projects/NeuralTF/scripts/visualize_results.py
```

Outputs 12 PNGs to `projects/NeuralTF/figures/`.

---

## 2. Building data from raw sources

These are prebuilt in the repo. Run only if upstream raw files change.

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
  --mmc7 "datasets/raw/Supplementary_Data_ King_2024/1-s2.0-S2211124724001712-mmc7.xlsx" \
  --out projects/NeuralTF/data/king_atlas.tsv
```

### 2.3 Raw data download reference

| Dataset | GEO Accession | Files needed |
|---------|--------------|--------------|
| Fincher 2018 | GSE111764 | `*_ClusteringDigitalExpressionMatrix.dge.txt.gz` (x3) |
| Plass 2018 | GSE103633 | `GSE103633_RAW.tar` + `*contigs.fasta.bz2` |
| King 2024 | Cell Reports | 8 supplementary xlsx: mmc2 through mmc8 |

The preprocessed h5ad files in `datasets/processed/` were built from the raw DGE files using one-time conversion scripts. The h5ad files are tracked in git for reproducibility.

---

## 2. Architecture

### Layer stack

```
Layer 0-2   Foundational: config, logging, exceptions, plugins
Layer 3     Sci stack: numpy, scipy, pandas
Layer 4     Bioinformatics: scanpy, anndata, leidenalg, igraph
Layer 5     CLI: click-based command interface
Layer 6     AI: OpenAI-compatible provider (stub by default)
Layer 7     Workflow: YAML-engine with step registry + provenance
Layer 8A    Omics: QC, normalize, cluster, trajectory, batch-correction
Layer 8B    Evidence: schema, scoring, cards, gene_mapping, bridge
Layer 8C    Ingest: auto-detect and load h5ad, DGE, TSV, 3X-mtx mixes
Layer 9     NeuralTF: pipeline.py + YAML workflow + demo-sim
Layer 10    UI: Streamlit app (Run / Results / Assistant)
```

### Evidence Integration Framework (8B)

```
EvidenceSource (Enum):
  EXPRESSION, SPECIFICITY, REPRODUCIBILITY,
  RNai, CORRELATION, FUNCTION,
  NEURAL_ENRICHED, NEURAL_SPECIFICITY

EvidenceRecord (dataclass):
  gene_id (e.g. dd_Smed_v6_15104_0_1)
  gene_name (e.g. soxB, coe)
  scores: dict[EvidenceSource -> float]
  notes: dict[EvidenceSource -> str]
  proof_status: str

EvidenceScorer:
  integrated_score(record) -> float  (renormalised weighted sum)

BridgeTable:
  v4_to_v6(), Process v6_to_v4(), v4|v6>= name()
```

### Pipeline call sequence

```
NeuralTFPipeline.run()
  load_datasets()          -> adata_fincher, adata_plass
  load_reference_tables()  -> tf_catalog, rnai_table, correlations, bridge
  run_qc()                 -> sweep QC + log1p + HVG + PCA + leiden clustering
  score_atlases()          -> per-cluster wilcoxon DE + expr(px_specificity)
  integrate_king_atlas()   -> seed neural-enriched TFs from G0 data, score expr/spec/neural
  integrate_rnai()         -> clean _embed matches via short dd-ID resolution
  integrate_correlations() -> match G0Â  GF pair gains against mmc6
  assign_reproducibility() -> n_atlases_supporting / 3
  write_outputs()          -> CSV + JSON + cards + AI in summary
```

---

## 3. Science behind NeuralTF

King et al. 2024 demonstrated that neural fate diversification in planarians happens overwhelmingly **post-mitotically**: ~7 broad neoblast clusters produce ~77 distinct G0 progenitor subclusters mapping onto 70+ known neuron types (glutamatergic, cholinergic, serotonergic, dopaminergic, GFI1B-, visual Anterior, etc.).

The prebuilt `king_atlas.tsv` captures which TF enriches for which G0 neural subcluster, with log2FC values from 1.5 to 8.0+. A TF with `log2FC >= 2.0` in a `neural_*` subcluster is considered **neural tissue.

We add two thesis-specific evidence streams:
- **NEURAL_ENRICHED (0.15 weight)**: 1.0 if the TF hits any `neural`* but is fully G0 subcluster with strong log2FC. 0.0 otherwise.
- **NEURAL_SPECIFICITY (0.10 weight)**: 1 / n_unique_neural_subclusters. A TF hitting 5 topics gets 0.20 (less specific). A TF precisely hitting 1 neural subtype gets 1.0 (best candidate for RNAi spam-stupidity).

This directly answers the thesis question: "Which TF, if knocked down by RNAi, would ablate a specific neuron subtype?"

### Evidence streams and weights

| Stream | Weight | Interpretation |
|--------|--------|----------------|
| Expression | 0.20 | max log2FC/5 across all 3 atlases |
| Specificity | 0.10 | 1 / n_clusters supporting the TF |
| Reproducibility | 0.15 | Fraction (n_atlases_supporting / 3) |
| RNAi | 0.15 | 1 if in King mmc5 table, 0 otherwise |
| Correlation | 0.10 | (G0_corr - X1_corr) × 3, capped to 1 |
| Function | 0.05 | 1 if known literature neural TF |
| Neural Enriched | 0.15 | 1 if neural G0 subcluster with log2FC >= 2 |
| Neural Specificity | 0.10 | 1 / n_unique_neural_subclusters |

### Tier assignment

- **HIGH**: RNAi-validated OR (streams ≥ 3 AND score ≥ 0.45)
- **MEDIUM**: streams ≥ 2 AND score ≥ 0.25
- **LOW**: everything else

### Proof status

| Status | Meaning | Priority |
|--------|---------|----------|
| `known_rnai_validated` | Already tested by King's own RNAi | Low — known |
| `novel_candidate` | Not in mmc5, no prior literature | **HIGH — run RNAi experiment** |
| `prior_fstf_not_tested` | Known FSTF, not yet RNAi found | Verified — re-test |

---

## 4. Common operations

### Find a gene by name across the bridge

```python
from bioforge.evidence import load_bridge
bridge = load_bridge('projects/NeuralTF/data/bridge.csv')
v6ids = bridge.df.loc[bridge.df['gene_name']=='soxB', 'v6_id']
print(v6ids.tolist())
```

### Query evidence data by e VI-id

```bash
python -c "
import pandas as pd
king = pd.read_csv('projects/NeuralTF/data/king_atlas.tsv', sep='\t')
hits = king[king['v6_id'] == 'dd_Smed_v6_10320_0_1']
print(hits[['subcluster','cell_type','log2fc']].to_string())
"
```

### Run a specific test

```bash
python -m pytest tests/unit/test_evidence.py::test_evidence_scorer_term\s -v
```

### Run all tests except optional deps

```bash
python -m pytest tests/unit/test_evidence.py tests/unit/test_evidence_cards.py tests/unit/test_evidence_readers.py -v
```

---

## 5. Directory map

| Path | What |
|------|------|
| `README.md` | Quick start |
| `bioforge.md` | This file — operations + architecture + science |
| `AGENTS.md` | AI-agents guide (subset of bioforge.md) |
| `pyproject.toml` | Package config, deps, CLI entry points |
| `src/bioforge/` | Python package root |
| `src/bioforge/evidence/` | Scoring, cards, gene mapping |
| `src/bioforge/projects/neuraltf/pipeline.py` | Main pipeline entry point |
| `projects/NeuralTF/data/` | bridge.csv, king_atlas.tsv |
| `projects/NeuralTF/scripts/visualize_results.pypy` | Figure generator |
| `projects/NeuralTF/runs/pipeline_run/` | Default output directory |
| `scripts/` | One-time build utilities |
| `tests/` | Test suite |
| `docs/` | Architecture decisions, development guide |
| `docker/` | Container setup |
| `datasets/processed/` | Preprocessed h5ad files | 

## 6. Reproducibility guarantees

- Explicit bridge table — no numeric-p Galileo ID guessing
- Deterministic random seed (42) per atlas subsampling
- Same code + same inputs = same outputs
- Locked dependencies in `docker/base/requirements.lock.txt`
- AI is optional — defaults to deterministic stub