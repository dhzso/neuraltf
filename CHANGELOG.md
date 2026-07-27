# Changelog

All notable changes to BioForge will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0-alpha.1] - 2026-07-22

### Added

- Repository initialization with foundational directory structure
- `README.md` with project overview, vision, goals, and repository layout
- `CHANGELOG.md` following the Keep a Changelog format
- `LICENSE` file under the MIT License
- `.gitignore` covering Python, Jupyter, IDE, operating system, Docker, and project-specific artifacts
- `docs/ARCHITECTURE.md` describing the planned layered architecture
- `docs/ROADMAP.md` describing development phases from Foundation to Open-source Release
- `docs/DEVELOPMENT.md` describing repository workflow, commit conventions, branch strategy, coding standards, versioning, and code review process
- `pyproject.toml` with PEP 621 project metadata, setuptools build backend, and src-layout configuration
- `src/bioforge/__init__.py` package initialization exposing `__version__`
- `tests/` directory structure with `unit/` and `integration/` subdirectories
- `docs/decisions/` directory for future Architecture Decision Records
- `.dockerignore` excluding `.git`, `datasets/`, `backups/`, Python bytecode, IDE, and Jupyter checkpoints from Docker build context
- `docker/base/Dockerfile` based on `python:3.11-slim-bookworm` with build toolchain, non-root user `bioforge` (UID/GID 1000), JupyterLab, and editable BioForge install
- `docker/scripts/` directory for future container lifecycle scripts
- `docker-compose.yml` with single `bioforge` service (`bioforge-dev` container), port 8888, repository bind mount, and `unless-stopped` restart policy
- ADR-0001 — Updated layer plan: split Layer 8 into 8A (Generic Omics Workflows) and 8B (Evidence Integration Framework), NeuralTF-directed
- Layer 3 Scientific Python Foundation with pinned versions: numpy 2.4.6, scipy 1.17.1, pandas 3.0.5, scikit-learn 1.9.0, matplotlib 3.11.1, plotly 6.9.0, statsmodels 0.14.6
- Layer 4 Bioinformatics Foundation added as `[bio]` optional extra: scanpy 1.11.5, anndata 0.12.6, scvelo 0.3.4, cellrank 2.0.7, biopython 1.87, gseapy 1.3.1, igraph 1.0.0, leidenalg 0.12.0, harmonypy 2.0.0
- `docker/base/requirements.lock.txt` — reproducible lock file capturing all 145 transitive dependencies at exact pinned versions
- Dockerfile updated to install with `pip install -e ".[bio]"` extra
- Layer 5 BioForge Core Framework:
  - `bioforge.core` — config (YAML loader + env overrides), logging, datasets (path resolution + traversal guard), exceptions hierarchy
  - `bioforge.plugins` — plugin base protocol + entry-point discovery + soft-fail loader
  - `bioforge.projects` — research project manager (NeuralTF-style scaffold validation + creation)
  - `bioforge.cli` — Click CLI with `info`, `datasets`, `projects`, `plugins` subcommands; `bioforge` console script registered
- 50 unit tests + 9 integration tests (against the real mounted repo) — all passing
- `click` and `pyyaml` added as runtime deps; `pytest` and `pytest-cov` added as `[dev]` extra

### Layer 8A — Generic Omics Workflows

- `bioforge.omics.qc` — `compute_qc_metrics`, `filter_cells`, `filter_genes`, `run_qc`; mt-prefix is optional (planarians have no MT chromosome)
- `bioforge.omics.normalize` — total-count normalize, log1p, highly-variable-gene flagging, `run_normalize` convenience pipeline
- `bioforge.omics.cluster` — PCA, neighbors graph, Leiden clustering, UMAP embedding
- `bioforge.omics.trajectory` — PAGA, scVelo velocity, CellRank 2.0 terminal-state estimation (CellRank 2.0.7 API: `compute_schur(n_components=...)`, `set_terminal_states`, explicit `compute_fate_probabilities`)
- `bioforge.omics.batch` — direct `harmonypy 2.0` batch-correction wrapper (avoids the scanpy wrapper's transpose bug where `Z_corr` is `(n_cells, n_comps)`)
- 12 unit tests covering all 8A modules; full suite 84 passed (50 core + 22 integration + 12 omics)
- `datasets/MANIFEST.md` documenting the untracked raw/reference files
- `.gitignore` updated to exclude `datasets/raw/` and `datasets/reference/`

### Layer 8B — Evidence Integration Framework (NeuralTF-directed)

- `bioforge.evidence.schema` — `EvidenceRecord`, `EvidenceSource`, `ConfidenceTier` dataclasses/enums
- `bioforge.evidence.gene_mapping` — `BridgeTable`, `load_bridge`, `build_bridge_from_names` (explicit v4↔v6 bridging; no numeric-prefix guessing)
- `bioforge.evidence.harmonization` — `AtlasHarmonizer` aligning Fincher/Plass/King cluster labels onto 9 canonical planarian tissue classes
- `bioforge.evidence.scoring` — `EvidenceScorer` weighted-sum of six evidence streams with renormalization for missing sources; `rank_candidates` returns sorted top-N
- `bioforge.evidence.confidence` — `assign_tiers` (high/medium/low) via `ConfidencePolicy` (≥4 streams ∧ ≥0.6 score for HIGH; ≥3 ∧ ≥0.35 for MEDIUM)
- `bioforge.evidence.ontology` — minimal `annotate_function` mapping known planarian TFs (soxB, pou4l-1, myoD, hnf4, …) to functional categories; stub interface for future GO mapping
- `bioforge.evidence.readers.king` — King 2024 mmc4 (TF catalog) / mmc5 (RNAi) / mmc6 (neural TF pair correlations) / mmc7 (G0 + X1 TF atlases) xlsx readers
- `bioforge.evidence.readers.fincher` & `bioforge.evidence.readers.plass` — load DGE matrices into `AnnData` (genes↔cells orientation corrected)
- ADR-0002 documents the module design, scoring formula, weight defaults (summing to 1.0), confidence thresholds, and independence from NeuralTF project specifics
- `openpyxl>=3.1,<4.0` added to `[bio]` extra
- 31 unit tests (18 evidence core + 13 readers); all hermetic (synthetic xlsx fixtures, no dependency on git-ignored raw datasets)
- Full test suite now **115 passing** (50 core + 9 integration + 12 omics + 31 evidence + 13 misc-refactor)
- `docker/base/requirements.lock.txt` refreshed with 154 pinned transitive deps (now includes openpyxl 3.1.5 and et_xmlfile 2.0.0)

### Layer 6 — AI Layer

- `bioforge.ai` — provider-agnostic AI interface (OpenAI-compatible HTTP)
- `bioforge.ai.OpenAICompatClient` posts to `${base_url}/chat/completions`
  (works with OpenAI, OpenRouter, Together, Groq, Anyscale, vLLM, SGLang,
  LM Studio, Ollama's OpenAI shim) via stdlib `urllib` (zero new runtime deps)
- `bioforge.ai.StubAssistant` deterministic fallback so workflows degrade
  gracefully when no `BIOFORGE_AI_API_KEY` is configured
- `bioforge.ai.tools` — `register_tool` decorator and built-in
  `lookup_gene` + `summarize_candidates` for LLM-assistant agency
- Config discovery: env vars (`BIOFORGE_AI_BASE_URL`, `BIOFORGE_AI_API_KEY`,
  `BIOFORGE_AI_MODEL`) → `~/.config/bioforge/ai.yaml` → safe defaults

### Layer 7 — Workflow Engine

- `bioforge.workflow.engine` — declarative YAML workflows with
  step references (`$step_id.output` or `{"$step": ..., "$output": ...}`),
  provenance capture, and progress callbacks for CLI / UI
- `bioforge.workflow.registry.StepRegistry` — `@register("name")` decorator
  with input/output contract tracking for assistant introspection
- `bioforge.workflow.steps` — builtins: `ingest`, `qc`, `normalize`,
  `cluster`, `trajectory`, `batch_correct`, `evidence.demo_rank`,
  `evidence.write_rank_csv`, `evidence.build_cards`, `report.write_cards_md`,
  `ai.summarize_candidates`
- New CLI: `bioforge run workflow.yaml --out runs/<ts>` executes workflows
  end-to-end; writes `provenance.json` + `summary.json` per run

### Layer 8C — Dataset Ingestion

- `bioforge.ingest` — `ingest_dataset(source)` accepts GEO/SRA accessions,
  URLs, or local paths and auto-detects format by **content sniffing**
  (HDF5 magic / gzip magic / text content), not extension alone
- Format families: `.h5ad`, `.csv[.gz]`, `.tsv[.gz]` / `.txt[.gz]` DGE,
  10x `matrix.mtx` directories
- Graceful `UnknownFormatError` so UI/CLI can render a friendly message
  instead of a stack trace (per ADR-0003); accessions defer to user-side
  geo/sra-toolkit fetch with a clear error message
- Optional FASTQ→matrix recipe declaration in
  `bioforge.ingest.fastq` (stub for future `[fastq]` extra)

### Layer 8B Extension — Cross-Atlas Evidence Cards

- `bioforge.evidence.cards` — `EvidenceCard`, `ProofStatus`,
  `build_evidence_card`, `build_cards_for_records`, `render_card_markdown`,
  `render_cards_markdown`
- Proof status separates **novel candidates** (no RNAi, prior unknowns —
  the thesis-essential class Deepanshu will validate) vs **known RNAi
  validated** vs **prior FSTF not tested** — all three contribute to the
  same confidence score, no evidence stream is discarded
- Markdown rendering produces per-candidate "evidence cards" suitable for
  both console reports and the Streamlit UI

### Layer 9 — NeuralTF

- `projects/NeuralTF/` directory with README, workflows, and runs/
- `projects/NeuralTF/workflows/demo_pipeline.yaml` — synthetic-data workflow
  that runs the full evidence → cards → AI-summary chain end-to-end
  (verified by integration test)
- `.gitignore` updated to exclude `runs/` and `projects/*/runs/`

### Layer 10 — Streamlit UI

- `bioforge.ui.app` — three pages (Run / Results / Assistant) using Streamlit
- Delegates to Layer 7's `WorkflowExecutor` for actual execution so CLI and
  UI behave identically
- AI panel degrades to `StubAssistant` when unconfigured
- `[streamlit]` optional extra added (streamlit + plotly); container now
  includes streamlit 1.60.0

### ADRs

- ADR-0003 — AI Layer + Workflow Engine + Dataset Ingestion + UI design
- ADR-0001 layer-plan table updated: Layers 6, 7, 8C, 9, 10 marked Completed

### Test totals

- **165 tests passing** (was 115)
- 21 new AI Layer tests, 21 new workflow-engine tests, 18 new ingest tests,
  7 new evidence-cards tests, 2 new UI smoke tests, 2 new NeuralTF
  integration tests
- `docker/base/requirements.lock.txt` now pins 172 transitive deps
  (+streamlit 1.60.0 + altair 6.2.2 + pyarrow 24.0.0 + httptools 0.8.0 …)
