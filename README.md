# BioForge · NeuralTF

A reproducible pipeline for **planarian neural-fate-specific transcription factor**
discovery. Integrates three peer-reviewed single-cell RNA-seq atlases
(Fincher 2018, Plass 2018, King 2024) plus the King 2024 RNAi phenotype
screen and neural TF-pair correlation data, scoring each candidate TF on
7 evidence streams and flagging priority targets for RNAi validation.

---

## Quick start

```bash
# 1. Set up Python 3.11+ environment
python -m venv .venv
.venv\Scripts\activate                       # Windows
source .venv/bin/activate                    # Linux/Mac

pip install -e ".[bio,streamlit]"

# 2. Place the raw downloads in datasets/raw/ (see "Datasets" below),
#    then build everything locally with ONE command:
python scripts/generate_all.py                # fincher/plass h5ads -> bridge/atlas
                                              # -> pipeline -> planmine -> shortlist -> figures
# (offline?  python scripts/generate_all.py --skip-planmine)

# 3. Run the Streamlit UI
bioforge ui                                  # Streamlit UI - http://localhost:8501
```

Pipeline outputs are written to `projects/NeuralTF/runs/pipeline_run/`:

| File | Content |
|------|---------|
| `rank_neural.csv`       | Neural-enriched candidates with proof_status |
| `rank.csv`              | All ~249 candidates |
| `evidence_cards.md`     | Per-candidate evidence summary |
| `pipeline_results.json` | Machine-readable top 50 |

### Prioritization outputs (RNAi shortlist)

Run after the pipeline (step 4 above). Written to `projects/NeuralTF/results/`:

| File | Content |
|------|---------|
| `top10_neural_tfs_prioritized.csv` | 5 Track A + 5 Track B shortlist |
| `candidate_summary_report.md`      | Per-candidate evidence + wet-lab suggestion |

Supporting data generated to `datasets/processed/`:
`planmine_annotations.parquet` (long-format GO/domain/BLAST rows) and
`planmine_transcripts.fasta` (transcript sequences, one record per candidate
that has one). Both are outputs of `scripts/query_planmine.py`.

---

## What it does

The pipeline seeds TF candidates from a King 2024 G0 atlas TF catalog,
scores them on 7 evidence streams, then filters for neural-fate specificity.

### Atlases

| Atlas | Year | Cells | Role |
|-------|------|-------|------|
| Fincher | 2018 | 50,562 | Whole-animal cell-type atlas (dd_Smed_v4) |
| Plass   | 2018 | 37,507 | Independent replication atlas (dd_Smed_v6) |
| King    | 2024 | G0 progenitors | Neural ground truth: enrichment across 77 neural G0-progenitor subclusters (of 175 subclusters in the processed atlas; 955 at the paper's full resolution) |

### Evidence streams (7)

| Stream | Weight | How it's computed |
|--------|--------|-------------------|
| Expression         | 0.211 | max log2FC/5 across all 3 atlases |
| Specificity        | 0.105 | 1 / n_clusters supporting the TF |
| Reproducibility    | 0.158 | atlases_supporting / 3 |
| RNAi               | 0.158 | 1 if gene is in King mmc5 RNAi table |
| Correlation        | 0.105 | G0-X1 pair correlation gain x 3 |
| Neural Enriched    | 0.158 | 1 if enrichment for any neural G0 subcluster with log2FC >= 2.0 |
| Neural Specificity | 0.105 | 1 / n_neural_subclusters present in |

Weights sum to 1.0; the scorer renormalizes over the streams present per
candidate. Note: `neural_enriched` is the cohort-defining neural gate
restated — within the neural subset it is constant 1.0 and only
discriminates neural vs non-neural candidates in the full 249 list.

### Tier assignment

- **HIGH**: RNAi-validated OR (streams >= 3 AND score >= 0.45)
- **MEDIUM**: streams >= 2 AND score >= 0.25
- **LOW**: all others

### Proof status

- `known_rnai_validated` — already tested by King et al. via RNAi
- `novel_candidate` — not yet tested — priority for new experiments
- `prior_fstf_not_tested` — known FSTFs from literature without RNAi data

---

## Dual-track prioritization (top-10 RNAi shortlist)

After the pipeline ranks ~99 neural-enriched candidates, two scripts turn that
rank list into a concrete wet-lab shortlist:

```bash
python scripts/query_planmine.py          # 1. Fetch annotations from PlanMine
python scripts/prioritize_neural_tfs.py   # 2. Score + split into tracks
```

**Step 1 — PlanMine annotation** (`scripts/query_planmine.py`) queries the
PlanMine InterMine API for every `dd_Smed_v6_*` candidate and stores GO terms,
protein domains (PFAM/InterPro), cross-species BLAST hits and the transcript
sequence in `datasets/processed/planmine_annotations.parquet` +
`planmine_transcripts.fasta`. Internal implementation:
`src/bioforge/projects/neuraltf/planmine.py`. Requires internet access to
`planmine.mpibpc.mpg.de` on the first run; the parquet/FASTA are cached
locally afterwards, so the prioritization step works offline.

**Step 2 — transparent dual-track scoring** (`scripts/prioritize_neural_tfs.py`):

- merges the mmc4 TF catalog (gene symbol, human ortholog, known-TF flag) and
  the PlanMine annotations onto the 99 candidates;
- maps each gene to its v4 alias via `projects/NeuralTF/data/bridge.csv`,
  flagging `unique` / `ambiguous` mappings (no numeric-ID guessing);
- applies small, additive, documented bonuses: TF domain `+0.05`, neural GO
  `+0.03`, TF GO `+0.02`, human ortholog `+0.02`, RNAi-validated `+0.02`
  (capped at 1.0, each category counted once);
- **Track A** = top 5 RNAi-validated candidates; **Track B** = top-5 novel
  candidates that pass a *tangible TF identity* filter (a DNA-binding-domain
  hit in PlanMine or an mmc4 TF flag — no hypothetical factors without a
  domain);
- appends cross-stage dynamics: Plass X1 neoblast mean vs G0 progenitor
  log2FC (requires `plass_v6.h5ad`; add `--skip-x1` to omit).

### Dirichlet-robust prioritization (sensitivity analysis)

The fixed-weight composite score uses a single weight vector
`W = [0.211, 0.105, 0.158, 0.158, 0.105, 0.158, 0.105]`. To test whether the
shortlist is robust to plausible weight perturbations, run the Dirichlet
analysis:

```bash
python projects/NeuralTF/scripts/dirichlet_prioritize.py             # 1000 draws, k=40
python projects/NeuralTF/scripts/dirichlet_visualize.py              # 5 Nature-style figures
```

This samples **1000 weight vectors** from a Dirichlet distribution centered
on `W` (concentration `k=40` ≈ 40 pseudo-observations; ~95% of weight mass
within ±0.1 of defaults). One weight vector per draw is applied to all 99
candidates; the median integrated score across draws is the "Dirichlet-robust"
score. The same Track A/B selection logic is then applied.

**Outputs** (in `projects/NeuralTF/results/` and `projects/NeuralTF/figures/`):
- `results/dirichlet_top10_prioritized.csv` — full top-10 (5 Track A + 5 Track B)
- `results/dirichlet_overall_top10.csv` — track-based shortlist
- `results/dirichlet_candidate_summary_report.md` — per-candidate evidence
- `figures/fig_dirichlet_*.png` — 5 publication-quality figures (track A/B, scatter, combined, score-shift)

**Key finding:** Track B is identical under both methods (top-5 stable). Track A
differs by one candidate (dd14712 replaces dd13343) — dd13343 relies heavily on
the RNAi stream (score=1.0); when the Dirichlet sampler down-weights that
stream, dd14712 (more balanced evidence) rises above it. Score shifts across
the top-10 are tiny (±0.006), confirming high robustness.

### Dirichlet-uniform (non-informative) prioritization

The centered Dirichlet tests "what if the defaults are approximately right?".
A second analysis tests "what does the data itself say?" with **no prior
preference** for any weighting (`alpha_i = 1` for all 7 streams):

```bash
python projects/NeuralTF/scripts/dirichlet_uniform.py                # 1000 draws, α=1
python projects/NeuralTF/scripts/dirichlet_uniform_viz.py            # 2 3-way comparison figures
python projects/NeuralTF/scripts/dirichlet_uniform_full_figures.py   # 5 Nature-style figures (mirrors centered)
python projects/NeuralTF/scripts/dirichlet_method_comparison.py      # 4 3-way method-comparison figures
```

**Outputs** (`results/` for tabular, `figures/` for PNGs, all gitignored):

| Script | Output files |
|---|---|
| `dirichlet_uniform.py` | `results/dirichlet_uniform_top10.csv`, `dirichlet_uniform_overall_top10.csv`, `dirichlet_uniform_full_rank.csv`, `dirichlet_uniform_summary.txt` |
| `dirichlet_uniform_viz.py` | `figures/fig_dirichlet_uniform_vs_centered.png`, `fig_dirichlet_3way_comparison.png` |
| `dirichlet_uniform_full_figures.py` | `figures/fig_dirichlet_uniform_trackA_top5.png`, `fig_dirichlet_uniform_trackB_top5.png`, `fig_dirichlet_uniform_scatter.png`, `fig_dirichlet_uniform_combined.png`, `fig_dirichlet_uniform_score_shift.png` |
| `dirichlet_method_comparison.py` | `figures/fig_dirichlet_score_density.png`, `fig_dirichlet_rank_correlation.png`, `fig_dirichlet_score_volatility.png`, `fig_dirichlet_method_summary.png` |

**Key finding:** 8/10 overlap with fixed-weight top-10. Track A: **dd31784**
(Homeobox, multi-subcluster) replaces dd13343. Track B: **dd33456** replaces
dd11930. Candidates unique to uniform Dirichlet are **fundamentally robust**:
they score high under ANY weighting, not just under the defaults.

**Note on scoring:** The 3-way comparison figure y-axis label is
"Base integrated score (before composite bonuses)" because the bars compare
the base scores from each method (fixed_weight_score, dirichlet_median_score,
uniform_median_score). The composite bonuses (TF domain, GO terms, RNAi
phenotype, etc.) are applied on top of each method's base score for the
final ranking — see `prioritize.py:_composite_score()`.

All Dirichlet scripts are integrated into `scripts/generate_all.py` as
steps 10–15 (run `python scripts/generate_all.py` after the main pipeline to
regenerate everything end-to-end).

### Filter breakdown (249 → 96 → 99)

The candidate-count numbers come from **stream-based filters applied upstream
of the integrated score**, not from score thresholds:

| Count | Filter | Source |
|---|---|---|
| **249** | `p ≤ 0.05` in ≥1 cluster's differential expression | `pipeline.py:286-290` |
| **96** | `neural_specificity.notna()` (≥1 King neural subcluster hit) | `pipeline.py:390-396` |
| **99** | `neural_specificity ∪ (rnai > 0)` — final neural-filtered set | `pipeline.py:573-574` |

The 99 count = 96 (King neural signal) ∪ 3 (RNAi-validated only: dd16955,
dd6626, dd12317). The 3 RNAi-only candidates have no King-atlas neural
signal but are kept because their biological validation is independent of
King 2024.

Scoring rules live in `src/bioforge/projects/neuraltf/prioritize.py`, which are unit-tested. The Streamlit UI shows the same tables under the
**Prioritization** page (`http://localhost:8501`).

---

## Datasets (no data files are bundled — you generate them locally)

The repository intentionally contains **no data files** (see `.gitignore`).
To reproduce every artifact:

1. Download the raw source files and place them in the locations below
2. Run:

   ```bash
   python scripts/generate_all.py            # everything, incl. PlanMine (network)
   # or offline: python scripts/generate_all.py --skip-planmine
   ```

   `generate_all.py` runs the nine steps below in dependency order, each
   gated on its inputs (missing inputs → skipped with a reason, never aborts).
   Steps whose outputs already exist are **skipped** — pass `--force` after a
   code change to regenerate everything:
   fincher h5ad → plass h5ad → bridge.csv → king_atlas.tsv → pipeline →
   PlanMine parquet (network) → prioritization shortlist + report → 13 main
   figures → 4 GO supplementary figures. Step-by-step commands:

**Fincher (GSE111764):**

1. Download from <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111764>
2. Extract `GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz` to
   `datasets/raw/GSE111764_GEO_Fincher_atlas/`
3. `python scripts/convert_fincher.py`

**Plass (GSE103633):**

1. Download the `GSE103633_RAW.tar` supplementary file from <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103633>
2. Place it anywhere under `datasets/raw/` (e.g. `datasets/raw/GSE103633_GEO_Plass_atlas/`); the script auto-locates the tar
3. `python scripts/consolidate_plass.py`

**King 2024 supplementary (Cell Reports):**

1. Download mmc2-mmc7.xlsx from the Cell Reports paper supplementary
2. Place under `datasets/raw/Supplementary_Data_ King_2024/`

**Rosetta Stone gene-ID bridge (v4↔v6):**

1. Download from <https://planosphere.stowers.org/pub/analysis/rosetta/smed_20140614.mapping.rosettastone.2020/smed_20140614.mapping.rosettastone.2020.txt>
2. Place at `datasets/raw/smed_20140614.mapping.rosettastone.2020.txt`
3. `python scripts/build_bridge.py`

**Gene Ontology (go.obo):**

1. Download the current release from the Gene Ontology downloads — either the
   mirror at <https://current.geneontology.org/ontology/go.obo> or the OBO
   Foundry URL <http://purl.obolibrary.org/obo/go.obo> (plain text, ~40 MB)
2. Place it at `datasets/raw/go.obo` (required by the GO supplementary
   figures, see `bioforge.md` §1.8)

**PlanMine annotations (network):** `python scripts/query_planmine.py`
(queries the PlanMine web API for the run's candidates; the docs can only
show results of the past local run that produced them).

The SRA records for the raw reads are Plass: SRP117156 (BioProject
PRJNA403817) and Fincher: SRP135258 (BioProject PRJNA438083).

---

## Project layout

```
Bioinformatics/
├── pyproject.toml                            Package config + dependencies
├── README.md                                 This file
├── bioforge.md                               Extended operations + architecture
├── .streamlit/config.toml                    Streamlit config (no email prompt)
│
├── src/bioforge/                             Python package
│   ├── evidence/                             8-stream scoring engine
│   │   ├── schema.py                         EvidenceRecord, EvidenceSource
│   │   ├── scoring.py                        Weighted score integration
│   │   ├── confidence.py                     Tier assignment (HIGH/MEDIUM/LOW)
│   │   ├── cards.py                           Per-candidate evidence cards
│   │   ├── gene_mapping.py                   v4<->v6 bridge table
│   │   └── readers/                          Atlas dataset readers
│   ├── projects/neuraltf/pipeline.py         Main pipeline (NeuralTFPipeline)
│   │   ├── planmine.py                       PlanMine InterMine REST client + classifiers
│   │   └── prioritize.py                     Dual-track scoring (pure, unit-tested)
│   ├── omics/                                scRNA-seq operations
│   ├── workflow/                             YAML workflow engine
│   ├── ai/                                   AI assistant layer
│   ├── cli/                                  Command-line interface
│   ├── ui/                                   Streamlit UI
│   └── core/                                 Config, datasets, logging, plugins
│
├── datasets/
│   ├── processed/                            Built h5ads + planmine parquet/fasta
│   │                                        (generated by scripts, gitignored)
│   ├── raw/                                  Raw downloads (gitignored)
│   └── reference/                            Reference tables (gitignored)
│
├── projects/NeuralTF/
│   ├── data/
│   │   ├── bridge.csv                        v4<->v6<->gene_name Rosetta stone
│   │   └── king_atlas.tsv                    Prebuilt G0 enrichment data
│   ├── scripts/visualize_results.py          Generate published figures
│   ├── figures/                              12 generated PNGs
│   ├── results/                              top10_neural_tfs_prioritized.csv + summary report
│   └── runs/                                 Pipeline output directory
│
├── scripts/                                   Utility scripts
│   ├── convert_fincher.py                     Build fincher_subsample.h5ad
│   ├── consolidate_plass.py                   Build plass_v6.h5ad
│   ├── build_bridge.py                        Rebuild bridge.csv
│   ├── build_king_atlas.py                    Rebuild king_atlas.tsv from mmc7
│   ├── audit_king_atlas.py                    Diagnostic: King atlas stats
│   ├── query_planmine.py                      Fetch PlanMine annotation/FASTA for candidates
│   ├── prioritize_neural_tfs.py               Build dual-track top-10 shortlist
│   └── run.py                                Pipeline runner (alternate entry)
│
├── tests/                                     Test suite (195 tests passing)
└── docs/                                      Architecture decisions
```

---

## CLI reference

```bash
bioforge --version
bioforge --help

# Main commands
bioforge neuraltf run [--subsample N] [--out DIR]   Run the NeuralTF pipeline
bioforge ui [--port 8501] [--host localhost]       Launch the Streamlit UI
    (pages: Run / Results / Prioritization / Assistant)

# Supporting commands (kept for advanced use)
bioforge info                                      Show build + config
bioforge datasets list [--category raw|processed]   List datasets
bioforge projects list                              List research projects
bioforge plugins list                               List plugins
bioforge run WORKFLOW.yaml                          Execute a YAML workflow
```

---

## Optional dependencies

The default `[bio]` extra installs the slim stack the pipeline actually
needs: scanpy, anndata, biopython, gseapy, igraph, leidenalg, openpyxl.

Heavy scverse extras (`harmonypy`, `scvelo`, `cellrank`) are *intentionally*
excluded because their native BLAS/CMake builds fail on a fresh Windows
install. The `bioforge.omics.batch` and `bioforge.omics.trajectory` modules
import cleanly without them; calling the wrapper functions raises a clear
`ImportError` if you forgot to install the underlying package.

Install them only if you actually use the corresponding wrapper:

```bash
pip install harmonypy    # bioforge.omics.batch.run_harmony
pip install scvelo       # bioforge.omics.trajectory.velocity
pip install cellrank     # bioforge.omics.trajectory.cellrank_terminal_states
```

---

## Reproducibility

- Atlases are independent experiments from peer-reviewed papers
- Gene IDs are bridged via a mandatory Rosetta Stone CSV
- Subsampling uses `random_state=42` everywhere
- AI operations use a deterministic `StubAssistant` unless an API key is configured
- All 195 unit tests pass on a clean install; 14 more skip when optional deps absent

## License

MIT — see [LICENSE](LICENSE)
