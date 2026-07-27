# Architecture

> High-level architecture for BioForge.
>
> This document describes the planned system architecture only.
> No implementation details are included here.

---

## Overview

BioForge is a containerized bioinformatics workstation running on Windows via
Docker Desktop. The system is built as a stack of layers, each layer building
on the one beneath it. This separation keeps concerns isolated and allows
individual layers to evolve without destabilizing the whole platform.

---

## Layered Architecture

```
┌─────────────────────────────────────────────┐
│  Windows                                     │  Host operating system
├─────────────────────────────────────────────┤
│  Docker Desktop                              │  Container runtime
├─────────────────────────────────────────────┤
│  BioForge Base                               │  Debian Slim foundation image
├─────────────────────────────────────────────┤
│  Scientific Python Stack                     │  Core scientific libraries
├─────────────────────────────────────────────┤
│  Bioinformatics Layer                        │  Analysis pipelines and tools
├─────────────────────────────────────────────┤
│  AI Layer                                    │  AI-assisted analysis and memory
├─────────────────────────────────────────────┤
│  CLI                                         │  Command-line interface
├─────────────────────────────────────────────┤
│  Plugins                                     │  Extensible modules
└─────────────────────────────────────────────┘
```

---

## Layer Descriptions

### Windows

The host operating system. BioForge does not run natively on Windows; all
computation occurs inside containers.

### Docker Desktop

Provides the container runtime on Windows. Manages image builds, container
lifecycle, volume mounts for datasets, and network isolation.

### BioForge Base

The foundation container image built on Debian Slim. Establishes system-level
dependencies, locale settings, and a stable Python 3.11 runtime. All higher
layers inherit from this base.

### Scientific Python Stack

Core scientific computing libraries (NumPy, SciPy, pandas, scikit-learn, etc.)
installed on top of the BioForge Base. This layer is shared across all
bioinformatics and AI workloads.

### Bioinformatics Layer

Domain-specific analysis pipelines: single-cell RNA-seq, bulk RNA-seq, spatial
transcriptomics, ATAC-seq, gene regulatory network inference, and report
generation. Each pipeline is modular and independently testable.

### AI Layer

AI-assisted analysis capabilities including prompt management, persistent
memory, and automated analysis assistance. This layer interacts with the
Bioinformatics Layer to augment — not replace — researcher-driven workflows.

### CLI

Command-line interface that exposes BioForge functionality to the user. Acts
as the primary entry point for running pipelines, managing projects, and
invoking AI-assisted features.

### Plugins

Extensible module system allowing users to add custom analysis pipelines,
data formats, or tool integrations without modifying core BioForge code.
Plugins are loaded at runtime by the CLI.
