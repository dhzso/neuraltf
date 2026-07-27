# ADR-0001: Updated BioForge Layer Plan — Layer 8 Split and NeuralTF-Directed Layer 8B

- **Status:** Accepted
- **Date:** 2026-07-26
- **Deciders:** Deepanshu (Owner), OpenCode (acting Architecture Lead)
- **Supersedes:** BioForge Engineering Master Plan v1.0 (Layer 8 section)

## Context

The BioForge Engineering Master Plan v1.0 defined Layer 8 as a single "Analysis
Pipelines" layer. After reviewing Deepanshu's MS thesis materials (Master
slides, King 2024, Fincher 2018, Plass 2018) and the supplied datasets, two
realities became clear:

1. **The thesis requires two distinct kinds of work** that have different
   reuse profiles and should not be conflated:
   - **Generic omics operations** (QC, normalization, clustering, trajectory,
     batch correction) — these are reusable Scanpy/scVelo/CellRank operations
     applicable to ANY scRNA-seq project.
   - **Evidence integration across independent cell atlases** — this is the
     novel scientific contribution: harmonizing gene identifiers across
     Fincher/Plass/King, scoring TFs by multi-criteria evidence, ranking
     candidates. This framework is specifically built for cross-atlas TF
     prioritization but is designed to be reusable for future studies.

2. **The charter's modularity principle** ("professional software product,
   not a collection of scripts") dictates that reusable infrastructure and
   research-specific logic live in separate layers.

3. **Thesis deadline pressure** argues against fully generic abstractions.
   However, BioForge's charter explicitly states "Never optimize for
   short-term convenience at the cost of architecture."

## Decision

**Split Layer 8 into two sub-layers:**

### Layer 8A — Generic Omics Workflows

Reusable single-cell analysis operations built on Scanpy/AnnData/scVelo/CellRank:
- Quality control
- Normalization
- Clustering
- Trajectory inference
- Batch correction

These are designed to work on any AnnData object, not just the NeuralTF datasets.

### Layer 8B — Evidence Integration Framework (NeuralTF-directed)

The novel framework for cross-atlas transcription factor prioritization. Built
with Deepanshu's thesis (NeuralTF) as the first consumer and primary use case,
BUT organized as a reusable framework with six components:

1. **Gene mapping** — Resolve v4↔v6 identifier bridges (Fincher uses dd_Smed_v4;
   King/Plass use dd_Smed_v6). No numeric-prefix guessing. Requires an explicit
   bridge table.
2. **Atlas harmonization** — Align cell-type annotations across independent
   atlases (Fincher/Plass/King cells → canonical tissue/lineage labels).
3. **Evidence scoring** — Multi-criteria score integrating:
   - expression enrichment (per-cluster log2FC, p-values from King S6)
   - cell-type specificity
   - cross-dataset reproducibility (Fincher + Plass independent confirmation)
   - functional annotation / regeneration relevance
   - RNAi phenotype support (King S4)
4. **Candidate ranking** — Rank TFs by integrated score; output top candidates.
5. **Confidence estimation** — Tier candidates (high/medium/low) based on
   number of supporting evidence streams.
6. **Ontology mapping** — Map candidates to functional categories and GO terms.

### Updated Layer Plan (Full)

| Layer | Name | Status |
|-------|------|--------|
| 0 | Repository Foundation | Completed |
| 1 | Core Python Project | Completed |
| 2 | Development Infrastructure | Completed |
| 3 | Scientific Python Foundation | Completed |
| 4 | Bioinformatics Foundation | Completed |
| 5 | BioForge Core Framework | Completed |
| 6 | AI Layer | Not started |
| 7 | Workflow Engine | Not started |
| 8A | Generic Omics Workflows | Completed |
| 8B | Evidence Integration Framework | Completed |
| 9 | Research Modules — NeuralTF | Not started |
| 10 | Reporting | Not started |
| 11 | Open-source Release | Not started |

## Alternatives Considered

1. **Single Layer 8 (no split)** — Rejected. Conflates generic infrastructure
   with the novel scientific contribution. Violates the single-responsibility
   principle and complicates testing.

2. **Fully generic Layer 8B** — Considered but rejected for now. Would delay
   thesis delivery by ~2-3 weeks for abstractions that have no immediate
   second consumer. The NeuralTF-directed approach is a pragmatic compromise:
   BioForge gets the reusable framework, Deepanshu gets the thesis results,
   and the framework is generalizable later (the bridge table, evidence schema,
   and scoring formula are dataset-agnostic).

3. **Embed 8B in NeuralTF (skip framework)** — Rejected. Tightly couples the
   scientific logic to the project directory. Impossible to reuse for other
   organisms/studies. Violates the charter's "professional software product"
   principle.

## Consequences

- **Positive:** Thesis-relevant work begins sooner. Framework is testable
  independently of NeuralTF. ADR establishes durability of the decision.
- **Negative:** Some NeuralTF-specific assumptions may leak into the framework.
  Mitigation: document assumptions in module docstrings; revisit at Layer 11
  (open-source release) for hardening.
- **Risk:** If a future consumer of Layer 8B has very different atlas structures,
  refactoring will be needed. Acceptable — the alternative (over-engineering
  now) is worse.

## Implementation Note

Layer 8B will live in `src/bioforge/evidence/` with submodules for each of
the six components. The NeuralTF project (Layer 9) lives in
`projects/NeuralTF/` and uses both 8A (for omics operations) and 8B
(for evidence integration) to produce the thesis results.
