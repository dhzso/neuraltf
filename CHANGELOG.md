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
