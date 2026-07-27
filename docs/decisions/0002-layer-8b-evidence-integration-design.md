# ADR-0002: Layer 8B Evidence Integration Framework — Module Design

- **Status:** Accepted
- **Date:** 2026-07-27
- **Deciders:** Deepanshu (Owner), OpenCode (acting Architecture Lead)
- **Relates to:** ADR-0001 (layer split)

## Context

ADR-0001 split Layer 8 into 8A (generic omics) and 8B (Evidence Integration
Framework, NeuralTF-directed). 8A is complete. This ADR records the concrete
module design for Layer 8B informed by inspecting the King 2024 supplementary
xlsx files (mmc2–mmc7) and the reference summary.

### What the King supplementary data actually contains

| File | Sheet(s) | Content |
|------|----------|---------|
| mmc2 | Sheet1 | Gene name → dd_Smed_v6 contig table (figure label reference) |
| mmc3 | Sheet1 | HMMER output: dd_Smed_v6 contigs vs. Pfam DBD models |
| mmc4 | All / TF / TF (additional) / BB / other | **TF catalog** — `Gene ID`, human BLAST hit, `TF?`, `FSTF?` flags |
| mmc5 | Sheet1 | **RNAi phenotype table** — `FSTF RNAi` × `cell-type markers tested`; col headers after row 3 are marker genes |
| mmc6 | Sheet1 | **Neural TF pair correlations** — `TF1`, `TF2`, `X1 Correlation`, `G0 Correlation`, `G0 Cluster` |
| mmc7 | 9 sheets | **TF Atlas**: G0 Progenitor TF Atlas / Pvalues / Log2FC, X1 TF Atlas / Pvalues / Log2FC, X1 Major Tissue Atlas / Pvalues / Log2F. Each "Atlas" sheet lists top-8 enriched FSTFs per subcluster with companion p-values and log2FC sheets. |

### Evidence streams for TF prioritization (NeuralTF thesis)

1. **Expression enrichment**: from mmc7 log2FC + p-values per subcluster
2. **Cell-type specificity**: from mmc7 atlas (which subclusters a TF appears in)
3. **Cross-dataset reproducibility**: re-compute enrichment in Fincher + Plass
   atlases (independent datasets) using identity-bridged genes — a TF that is
   neural-enriched in King AND Fincher AND Plass is strong.
4. **Functional / regeneration relevance**: human homolog annotation from mmc4.
5. **RNAi phenotype support**: from mmc5 (does RNAi of this TF cause a marker-loss phenotype?).
6. **Pairwise TF correlation**: from mmc6 (coordinated switching G0 vs X1 — supports combinatorial code).

### Gene-identifier bridge problem (central engineering concern)

- Fincher = **dd_Smed_v4**
- Plass + King = **dd_Smed_v6**
- No numeric-prefix guessing. We need an explicit bridge table.
- mmc2 provides gene-name ↔ dd_Smed_v6 contig mapping. Fincher's own
  supplementary (in `datasets/raw/GSE111764.../family.soft.gz`) lists v4
  contigs with names. We will build the bridge by **gene-name matching**
  (where names are unique) into a `BridgeTable` that the framework loads.

## Decision

### Module layout

```
src/bioforge/evidence/
  __init__.py
  schema.py        # EvidenceRecord dataclass; ConfidenceTier enum; Enums for evidence sources
  gene_mapping.py  # BridgeTable dataclass; load/save; build_bridge_from_names
  harmonization.py # AtlasHarmonizer; align cell-type labels to a canonical tissue vocabulary
  scoring.py       # multi-criterion EvidenceScorer; weighted sum + per-source normalization
  ranking.py       # rank_candidates; output top-N with integrated scores
  confidence.py   # tier assignment: high/medium/low based on number of supporting streams
  readers/
    __init__.py
    king.py        # load King mmc4 (TF catalog), mmc5 (RNAi), mmc6 (correlations), mmc7 (atlas)
    fincher.py     # load Fincher cluster DGE / cluster labels
    plass.py       # load Plass cluster labels and PAGA tree
  ontology.py     # (stubbed for now) future mapping to GO / functional categories
```

### Public API (summary)

```python
from bioforge.evidence import (
    EvidenceRecord, ConfidenceTier, EvidenceSource,
    BridgeTable, AtlasHarmonizer, EvidenceScorer, rank_candidates,
)
from bioforge.evidence.readers import king, fincher, plass
```

### Scoring formula

Each TF receives a normalized score in `[0, 1]` per evidence source, then a
weighted integrated score. Weights are configurable with sensible defaults:

| Source | Default weight | Normalization |
|--------|---------------|----------------|
| Expression enrichment (King mmc7 log2FC, p<0.001) | 0.30 | max-normalized log2FC |
| Cell-type specificity (1 – Shannon entropy over clusters) | 0.20 | min-max scaled |
| Cross-dataset reproducibility (Fincher ∧ Plass ∧ King) | 0.25 | binary recoded as 0/1 per dataset (max 3) |
| RNAi phenotype support (mmc5 marker loss) | 0.15 | binary |
| Pairwise TF correlation gain G0−X1 (mmc6) | 0.10 | z-scored |
| Functional annotation (human homolog) | 0.05 | binary (homolog present) |

Integrity rule: the integrated score is `Σ wᵢ · sᵢ` with `Σ wᵢ = 1`.

### Confidence tiers

- **high**: ≥ 4 supporting streams AND integrated score ≥ 0.6
- **medium**: 3 supporting streams AND integrated score ≥ 0.35
- **low**: anything else

### Bridge table format

A `BridgeTable` is a `pandas.DataFrame` with columns
`["gene_name", "v6_id", "v4_id"]`. Some v4_ids may be `NaN` if a v6 gene has
no Fincher counterpart. The bridge is loaded from a CSV path (`bioforge.evidence.gene_mapping.load_bridge`).

For NeuralTF, the bridge CSV will be produced by the project, not bundled
with the framework — but the framework ships a `build_bridge_from_names`
helper to construct it from a name-matching pass.

### Independence from NeuralTF specifics

The framework takes generic inputs (AnnData objects with `obs['tissue']`,
`obs['cluster']`; TF catalog as a DataFrame; RNAi table as a DataFrame). It
does not import any NeuralTF project code, keeping Layer 9 separable from
Layer 8B as required by ADR-0001.

## Alternatives Considered

1. **Single module** — rejected; six logically distinct components deserve
   separate modules for testability.

2. **Hardcode weights for NeuralTF** — rejected; weights must be
   overridable so the framework is reusable.

3. **Makeontology.py wait for Layer 9** — rejected; the framework defines
   the stub interface now so callers can extend it later without touching 8B
   internals.

## Consequences

- Adds runtime dependency on `openpyxl` (for reading King xlsx) — to be added
  to `pyproject.toml` `[bio]` extra.
- Adds unit tests per module (≈ 30 tests expected), bringing the full suite
  to ~115 tests.
- NeuralTF (Layer 9) will wire concrete dataset paths and the bridge CSV
  through this framework; the framework itself runs on AnnData/DataFrame inputs.
