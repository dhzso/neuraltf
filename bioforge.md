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

The processed h5ad files are not committed (they sum to ~80 MB). Build them
locally from raw GEO downloads:

```bash
python scripts/convert_fincher.py      # Fincher atlas (GSE111764)
python scripts/consolidate_plass.py    # Plass atlas (GSE103633)
```

The King 2024 supplementary xlsx files (mmc4-mmc7) must be placed under
`datasets/raw/Supplementary_Data_ King_2024/`. The pipeline auto-discovers
them — the exact Elsevier filename (`1-s2.0-S2211124724001712-mmcN.xlsx`) is
tried first, then any `*mmcN.xlsx` in the directory.

### 1.4 Run the NeuralTF pipeline (CLI)

```bash
bioforge neuraltf run [--subsample 0] [--out projects/NeuralTF/runs/my_run]
```

This runs the full pipeline and writes results:

| File | Content |
|------|---------|
| `rank.csv` | All ~249 TF candidates ranked by integrated score |
| `rank_neural.csv` | Neural-enriched candidate subset with proof_status |
| `evidence_cards.md` | Per-candidate markdown evidence summary |
| `pipeline_results.json` | Machine-readable top 50 with tier metadata |

### 1.5 Prioritize RNAi targets (post-run: PlanMine annotation + dual-track top-10)

After a successful run, two scripts turn the ranked candidates into a wet-lab
shortlist:

```bash
python scripts/query_planmine.py             # 1. PlanMine annotations (needs internet)
python scripts/prioritize_neural_tfs.py      # 2. dual-track top-10 -> results/
```

**Step 1 — `scripts/query_planmine.py`**: queries the PlanMine InterMine API
(`https://planmine.mpibpc.mpg.de/planmine/service`) for every `dd_Smed_v6_*`
candidate in `rank_neural.csv` and writes:

| File | Content |
|------|---------|
| `datasets/processed/planmine_annotations.parquet` | Long-format rows: GO terms, protein domains (PFAM/InterPro), cross-species BLAST hits |
| `datasets/processed/planmine_transcripts.fasta` | Transcript sequence per candidate (dsRNA/FISH design) |

The client (`src/bioforge/projects/neuraltf/planmine.py`) handles retries,
rate limiting (~0.25 s/request) and logs 400 responses. `dd_Smed_v6_*` IDs are
*Contig* (Transcript subclass) records in PlanMine, so one query per gene per
annotation type is used — unauthenticated InterMine IN-constraints fail past
~400 items, per-contig queries stay well under that. Both outputs are kept
locally (gitignored), so re-runs of command 2 run offline against them; only
the first time after a fresh clone needs the network.
Smoke-test with `--limit 5`; use `--out`/`--fasta` to redirect outputs.

**Step 2 — `scripts/prioritize_neural_tfs.py`**: pure scoring logic in
`src/bioforge/projects/neuraltf/prioritize.py`:

- merges the mmc4 TF catalog (name/ortholog/TF flag) with PlanMine annotations;
- maps v4 aliases via `projects/NeuralTF/data/bridge.csv`, tagging
  `unique` / `ambiguous` (multi-hit rows get a blank v4 ID — no guessing) /
  `unmapped`;
- additive composite bonuses, each category counted once, capped at 1.0:
  DNA-binding TF domain `+0.05`, neural GO `+0.03`, TF GO `+0.02`, human
  ortholog `+0.02`, RNAi-validated `+0.02`;
- **Track A** = top-5 RNAi-validated targets; **Track B** = top-5 novel
  targets passing a TF-identity gate (PlanMine DNA-binding-domain hit or mmc4 TF flag);
- cross-stage dynamics: Plass X1 neoblast mean vs G0 progenitor log2FC
  (needs `plass_v6.h5ad`; add `--skip-x1` to omit).

Outputs (generated, gitignored):

| File | Content |
|------|---------|
| `projects/NeuralTF/results/top10_neural_tfs_prioritized.csv` | 5 Track A + 5 Track B rows incl. all bonus/dynamic columns |
| `projects/NeuralTF/results/candidate_summary_report.md` | Per-candidate evidence + wet-lab suggestion (dsRNA/FISH nt range) |

Only step 1 needs internet. The `mmc4`/`mmc5` xlsx under
`datasets/raw/Supplementary_Data_ King_2024/` are auto-discovered (see 2.2).

### 1.5b Dirichlet-robust prioritization (weight sensitivity)

The fixed-weight composite uses one weight vector
`W = [0.211, 0.105, 0.158, 0.158, 0.105, 0.158, 0.105]`. To test whether the
shortlist is robust to plausible weight perturbations, run the Dirichlet
analysis:

```bash
python projects/NeuralTF/scripts/dirichlet_prioritize.py
python projects/NeuralTF/scripts/dirichlet_visualize.py
```

Or just run the full one-command pipeline (`scripts/generate_all.py`) — the
four Dirichlet steps are integrated as steps 10–13.

**Method:** Samples 1000 weight vectors from `Dirichlet(alpha = W × k)` with
`k=40` (≈ 40 pseudo-observations; ~95% of weight mass within ±0.1 of
defaults). One weight vector per draw is applied to **all 99 candidates**;
NaN streams are zeroed. The **median** integrated score across draws is the
"Dirichlet-robust" score. Same Track A/B selection logic.

**Outputs** (`results/` for CSVs/MD, `figures/` for PNGs, all gitignored):
- `results/dirichlet_top10_prioritized.csv` — full top-10 (5 Track A + 5 Track B)
- `results/dirichlet_overall_top10.csv` — track-based shortlist
- `results/dirichlet_overall_top10_byscore.csv` — overall top-10 by score
- `results/dirichlet_candidate_summary_report.md` — per-candidate evidence + comparison table
- `figures/fig_dirichlet_*.png` — 5 publication-quality figures (Nature style):
  - `fig_dirichlet_trackA_top5.png` — Track A horizontal bars
  - `fig_dirichlet_trackB_top5.png` — Track B horizontal bars
  - `fig_dirichlet_scatter.png` — robustness: fixed vs Dirichlet
  - `fig_dirichlet_combined.png` — both tracks with composite bonus overlay
  - `fig_dirichlet_score_shift.png` — Dirichlet − fixed per candidate

**Result:** Track B identical under both methods. Track A: dd14712 (Ets,
balanced evidence) replaces dd13343 (Homeobox, RNAi-dependent) — the only
swap. Score shifts ±0.006 confirm high robustness.

### 1.5c Dirichlet-uniform (non-informative) prioritization

The centered Dirichlet above tests "what if the defaults are approximately
right?". The uniform Dirichlet tests "what does the data itself say?" with
**no prior preference** for any weighting:

```bash
python projects/NeuralTF/scripts/dirichlet_uniform.py
python projects/NeuralTF/scripts/dirichlet_uniform_viz.py
python projects/NeuralTF/scripts/dirichlet_uniform_full_figures.py
python projects/NeuralTF/scripts/dirichlet_method_comparison.py
```

**Method:** Uses `alpha_i = 1` for all 7 streams (uniform over the 7-simplex).
Samples 1000 weight vectors uniformly; one vector per draw applied to all 99
candidates; median score is the "uniform-robust" score.

**Outputs** (gitignored):
- `results/dirichlet_uniform_top10.csv` — track-based 5A+5B under uniform prior
- `results/dirichlet_uniform_overall_top10.csv` — overall top-10 by uniform median
- `results/dirichlet_uniform_full_rank.csv` — all 99 candidates with both scores
- `results/dirichlet_uniform_summary.txt` — 3-way comparison stats
- `figures/fig_dirichlet_uniform_vs_centered.png` — scatter (centered vs uniform)
- `figures/fig_dirichlet_3way_comparison.png` — grouped bars (fixed/centered/uniform)
  - **Y-axis label**: "Base integrated score (before composite bonuses)" — the
    bars compare base scores from each method; composite bonuses are applied
    separately for the final ranking (see `prioritize.py:_composite_score()`)
- `figures/fig_dirichlet_uniform_trackA_top5.png` — Track A (uniform)
- `figures/fig_dirichlet_uniform_trackB_top5.png` — Track B (uniform)
- `figures/fig_dirichlet_uniform_scatter.png` — fixed vs uniform
- `figures/fig_dirichlet_uniform_combined.png` — both tracks with composite
- `figures/fig_dirichlet_uniform_score_shift.png` — uniform − fixed per candidate
- `figures/fig_dirichlet_score_density.png` — KDE of all 99 for 3 methods
- `figures/fig_dirichlet_rank_correlation.png` — Spearman ρ heatmap
- `figures/fig_dirichlet_score_volatility.png` — per-candidate score range
- `figures/fig_dirichlet_method_summary.png` — 4-panel summary (overlap, tracks, volatility, composite effect)

**Result:** 8/10 overlap with fixed-weight. Track A: **dd31784** (Homeobox,
multi-subcluster evidence) replaces dd13343. Track B: **dd33456** replaces
dd11930. Candidates unique to uniform Dirichlet are **fundamentally robust**
— they score high under ANY weighting.

### 1.5d Filter breakdown (249 → 96 → 99)

The candidate-count numbers in this pipeline come from **stream-based
filters applied upstream of the integrated score** — not from score
thresholds. Here is the exact funnel:

| Step | Count | Filter | Where |
|---|---|---|---|
| All TF targets | ~2,800 | PlanMine + mmc6 (gene is a TF) | `pipeline.py:302-307` |
| **Expression p-value** | **249** | `best_p ≤ 0.05` in ≥1 cluster's differential expression | `pipeline.py:286-290` |
| Reproducibility + Specificity | 249 | streams populated automatically | `pipeline.py:380-385, 504-512` |
| RNAi stream | 249 | `rnai > 0` from King 2024 mmc5 | inherited from Fincher/Plass |
| King neural enrichment | 96 | `neural_specificity.notna()` (≥1 hit in any of King 2024's 77 neural subclusters) | `pipeline.py:390-396` |
| **Final neural-filtered set** | **99** | `neural_specificity.notna() ∪ (rnai > 0)` | `pipeline.py:573-574` |

**Key point:** the 99 count is `neural_specificity ∪ rnai`, not pure
`neural_specificity`. The 3 candidates that come from the `rnai > 0` branch
alone (dd16955, dd6626, dd12317) are **known RNAi-validated TFs** with high
integrated scores but no King-atlas neural signal — they are included
because their biological validation is independent of King 2024.

### 1.6 Launch the Streamlit UI

```bash
bioforge ui [--port 8501] [--host localhost]
```

Opens http://localhost:8501. Pages: **Run** (dataset status + live pipeline
run), **Results** (rank CSVs + full Visualization panel with 13 figures),
**Prioritization** (coverage metrics + Track A/B tables + the full markdown
report), **Assistant** (AI chatbot, see next section).

### 1.7 Generate publication figures (as PNG)

```bash
python projects/NeuralTF/scripts/visualize_results.py
```

Outputs 12 PNGs to `projects/NeuralTF/figures/`.

### 1.8 Generate the GO supplementary figures

```bash
python projects/NeuralTF/scripts/make_supp_go_figures.py
```

Outputs the 4 supplementary GO figures + the reduced gene×term matrix CSV to
`projects/NeuralTF/figures/supplementary/`. The script resolves every PlanMine
GO term against `datasets/raw/go.obo` (canonical names/namespaces, obsolete
terms excluded), so the figures always agree with the +0.03 neural-GO composite
bonus in the prioritization report. Pass `--help` to override the run dir,
parquet, top-10 shortlist or go.obo path.

**Prerequisite — go.obo:** download the current release to
`datasets/raw/go.obo` (a fresh clone has no bundled data): the Gene Ontology
downloads (<https://current.geneontology.org/ontology/go.obo> or
<http://purl.obolibrary.org/obo/go.obo>, ~40 MB plain text). The script
warns and falls back to PlanMine's term names if the file is missing. Full
context: `datasets/MANIFEST.md` → "Gene Ontology (go.obo)".

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
| Expression | 0.211 | max log2FC/5 across all 3 atlases |
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
