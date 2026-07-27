# BioForge

An open-source AI-native bioinformatics workstation.

---

## Project Overview

BioForge is a reproducible computational biology platform. The initial objective is to support single-cell RNA sequencing, transcription factor discovery, and planarian regeneration research.

The long-term objective is to evolve into a modular computational biology platform supporting:

- Single-cell RNA-seq
- Bulk RNA-seq
- Spatial transcriptomics
- ATAC-seq
- Gene regulatory network inference
- Automated AI-assisted analysis
- Report generation
- Plugin-based extensions

---

## Vision

Treat BioForge as a professional software product rather than a collection of scripts. Every file has a clear purpose. Every architectural decision supports long-term maintainability and reproducible scientific research.

---

## Goals

- **Stability** — predictable, reliable behavior across environments
- **Reproducibility** — identical results from identical inputs
- **Documentation** — every component is explained and discoverable
- **Modularity** — components are independently testable and replaceable
- **Maintainability** — the codebase is readable and extensible

Never optimize for short-term convenience at the cost of architecture.

---

## Repository Overview

```
D:\Bioinformatics
│
├── ai/                    AI prompts, memory, and configurations
│   ├── prompts/
│   ├── memory/
│   └── configs/
│
├── backups/               Project backups (gitignored)
│
├── config/                Global configuration files
│
├── datasets/              Data storage
│   ├── cache/             Transient cache (gitignored)
│   ├── processed/         Cleaned and processed data
│   ├── raw/               Raw sequencing data
│   └── reference/         Reference genomes and annotations
│
├── docker/                Container definitions
│   ├── base/              Base image
│   ├── compose/           Docker Compose files
│   └── services/          Service-specific images
│
├── docs/                  Project documentation
│
├── projects/              Individual research projects
│   └── NeuralTF/          Neural transcription factor project
│       ├── data/
│       ├── docs/
│       ├── figures/
│       ├── logs/
│       ├── notebooks/
│       ├── results/
│       └── scripts/
│
├── src/                   BioForge library source code
│
├── templates/             Project and analysis templates
│
├── tools/                 Standalone utility scripts
│
└── workspace/             Active working area
```

---

## Planned Technology Stack

| Component            | Technology                          |
|----------------------|-------------------------------------|
| Programming language | Python 3.11                         |
| Base OS image        | Debian Slim                         |
| Containerization     | Docker / Docker Compose             |
| Version control      | Git                                 |
| Documentation        | Markdown                            |
| Code style           | PEP8 with type hints                |

---

## Installation

> **Placeholder** — Installation instructions will be provided in a future milestone once the containerized environment is established.

---

## Roadmap

> **Placeholder** — See [docs/ROADMAP.md](docs/ROADMAP.md) for the full development roadmap.

---

## License

MIT License — see [LICENSE](LICENSE).
