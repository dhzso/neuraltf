# BioForge · NeuralTF

A reproducible pipeline for **planarian neural-fate-specific transcription factor** discovery. Integrates **five single-cell and regulatory atlases** (Fincher 2018, Plass 2018, Cui 2023, King 2024, Perez 2025) with 8 evidence streams, Bayesian Dirichlet uncertainty quantification, and ANANSE gene regulatory network validation to prioritize high-confidence targets for RNAi and functional validation.

> **Latest Pipeline Status:** 289 candidate TFs scored across 8 evidence streams → 102 neural-enriched candidates → top-10 consensus across fixed-weight, centered Dirichlet, and uniform Dirichlet methods.

---

## Quick Start

```bash
# 1. Set up Python 3.12+ environment
git clone https://github.com/dhzso/neuraltf.git
cd neuraltf

python -m venv .venv
.venv\Scripts\activate                       # Windows
# source .venv/bin/activate                  # Linux/Mac

pip install -e ".[bio,streamlit]"

# 2. Place raw downloads in datasets/raw/ (see datasets/MANIFEST.md)
#    and build all processed artifacts with ONE master orchestrator:
python scripts/generate_all.py

# 3. Or run the pipeline and downstream analyses step-by-step:
python scripts/run.py                        # Core pipeline (Fincher, Plass, Cui, King, Perez)
python scripts/run_downstream.py             # Dirichlet UQ, ANANSE scan, tables & 21 figures

# 4. Launch the interactive Streamlit UI
bioforge ui                                  # http://localhost:8501
```

### Pipeline Outputs (`projects/NeuralTF/runs/pipeline_run/`)

| File | Content |
|------|---------|
| `rank.csv` | All **289 candidates** with scores across all 8 evidence streams |
| `rank_neural.csv` | **102 neural-enriched candidates** with proof status |
| `evidence_cards.md` | Per-candidate markdown evidence summary (289 cards) |
| `pipeline_results.json` | Machine-readable candidate metadata |
| `checkpoint_0*.parquet` | 6 incremental audit checkpoints (atlas loads, post-QC, post-scoring, King, Perez, stream matrix) |

### Prioritization & Analysis Outputs (`projects/NeuralTF/results/`)

| File | Content |
|------|---------|
| `master_tf_catalog.csv` | Unified King + Perez TF catalog (14,682 unique v6 TFs) |
| `dirichlet_centered_all249_full_rank.csv` | Full Dirichlet-centered (k=40) rank list (289 candidates) |
| `dirichlet_uniform_all249_full_rank.csv` | Full Dirichlet-uniform (α=1) rank list (289 candidates) |
| `dirichlet_top10_prioritized.csv` | Dual-track top-10 shortlist (5 Track A + 5 Track B) |
| `ananse_network_full.csv` | ANANSE GRN scan across all 289 candidates & 9 cell fates |
| `ananse_top_regulators.csv` | Top planarian neural regulators ranked by out-degree |
| `supplementary_table_S1` – `S5` | Method comparison, full rank tables, and Dirichlet analyses |

---

## What It Does

The pipeline seeds candidate TFs across five single-cell and regulatory atlases, computes an 8-stream multi-evidence matrix, evaluates scoring stability under Dirichlet uncertainty sampling, and maps candidates into cell-fate regulatory circuits.

### Atlases Integrated (5)

| Atlas | Year | Modality / Data Type | Role in Pipeline |
|-------|------|----------------------|------------------|
| **Fincher** | 2018 | scRNA-seq (50.5K cells, v4 IDs) | Whole-animal cell-type expression & neural specificity |
| **Plass** | 2018 | scRNA-seq (37.5K cells, v6 IDs) | Independent whole-anatomy replication & X1 dynamics |
| **Cui** | 2023 | scRNA-seq (55.0K cells, 8 time points) | High-resolution regeneration time-course expression |
| **King** | 2024 | Single-cell TF catalog + RNAi screen | G0/X1 cluster enrichment, RNAi phenotypes (mmc5), TF pair correlations (mmc6) |
| **Perez** | 2025 | Lineage atlas + ANANSE GRNs | Lineage TF classification (MOESM5) & ANANSE GRN target validation (MOESM22) |

---

## Evidence Streams & Scoring Model (8 Streams)

Scoring utilizes a transparent, weighted multi-evidence integration model. Weights renormalize over streams present for each candidate:

$$\text{Integrated Score} = \sum_{i \in \text{Present}} w_i \cdot s_i \Bigg/ \sum_{i \in \text{Present}} w_i$$

| # | Stream | Default Weight ($w_i$) | Biological Basis & Computation |
|---|--------|------------------------|---------------------------------|
| 1 | **Expression** | 0.200 | $\min(1.0, \max(\text{log}_2\text{FC})/5)$ across Fincher, Plass, and Cui scRNA-seq atlases |
| 2 | **Specificity** | 0.100 | $1 / n_{\text{clusters}}$ supporting differential expression |
| 3 | **Reproducibility** | 0.100 | $n_{\text{atlases supporting}} / 5$ (Fincher, Plass, Cui, King, Perez) |
| 4 | **RNAi** | 0.100 | Binary indicator (1.0) if functional phenotype observed in King mmc5 screen |
| 5 | **Correlation** | 0.100 | $\min(1.0, \Delta r_{\text{G0-X1}} \times 3.0)$ co-expression gain from King mmc6 |
| 6 | **Neural Enriched** | 0.100 | Binary indicator (1.0) for G0 neural subcluster log₂FC ≥ 2.0 (King mmc7) |
| 7 | **Neural Specificity** | 0.100 | $1 / n_{\text{neural subclusters}}$ present in King atlas |
| 8 | **Perez Lineage** | 0.100 | Perez 2025 lineage TF class: **1.0** for neural-class (bHLH, Homeobox, POU, C2H2, ETS, etc.), **0.5** for other TF classes, **0.0** if absent |

### Confidence Tiers & Proof Status

- **Tier Assignment**:
  - **HIGH**: RNAi-validated OR (supporting streams $\ge 3$ AND score $\ge 0.45$)
  - **MEDIUM**: supporting streams $\ge 2$ AND score $\ge 0.25$
  - **LOW**: All other candidates
- **Proof Status**:
  - `known_rnai_validated` — Confirmed neural/regeneration phenotype in King et al. RNAi screen
  - `novel_candidate` — Uncharacterized TF with strong multi-atlas support (priority for wet-lab knockout)
  - `prior_fstf_not_tested` — Documented fate-specifying TF (FSTF) from literature without RNAi data

---

## Top Prioritized Candidates (Current Run)

| Rank | Gene ID (v6) | Gene Name | Integrated Score | Dirichlet Median | Proof Status | Primary DNA-Binding Domain / Family |
|:---:|:---|:---|:---:|:---:|:---|:---|
| 1 | `dd_Smed_v6_2946_0_1` | **dd2946** | 0.852 | 0.857 | `known_rnai_validated` | C2H2 Zinc Finger |
| 2 | `dd_Smed_v6_16472_0_1` | **dd16472** | 0.790 | 0.792 | `known_rnai_validated` | Homeobox |
| 3 | `dd_Smed_v6_38342_0_1` | **dd38342** | 0.782 | 0.783 | `known_rnai_validated` | POU / Homeobox |
| 4 | `dd_Smed_v6_4048_0_1` | **dd4048** | 0.780 | 0.791 | `novel_candidate` | bHLH |
| 5 | `dd_Smed_v6_14115_0_1` | **dd14115** | 0.764 | 0.762 | `known_rnai_validated` | Homeobox / LIM |
| 6 | `dd_Smed_v6_11150_0_1` | **dd11150** | 0.760 | 0.764 | `known_rnai_validated` | C2H2 Zinc Finger |
| 7 | `dd_Smed_v6_14824_0_1` | **dd14824** | 0.756 | 0.762 | `known_rnai_validated` | C2H2 Zinc Finger |
| 8 | `dd_Smed_v6_19890_0_1` | **dd19890** | 0.753 | 0.753 | `known_rnai_validated` | T-box |
| 9 | `dd_Smed_v6_31217_0_1` | **dd31217** | 0.751 | 0.760 | `novel_candidate` | bHLH |
| 10 | `dd_Smed_v6_6626_0_1` | **dd6626** | 0.750 | 0.749 | `known_rnai_validated` | Nuclear Hormone Receptor |

> **Key Finding:** **10/10 overlap** in top-10 candidates across fixed-weight, centered Dirichlet, and uniform Dirichlet methods confirms that candidate ranking is highly robust to parameter choices.

---

## Uncertainty Quantification & Sensitivity (Dirichlet Sampling)

To test sensitivity against arbitrary weighting assumptions, we employ Monte Carlo Dirichlet weight perturbation (1,000 draws, seed=2024):

1. **Centered Dirichlet ($k = 40$)**:
   $$\mathbf{w}^{(m)} \sim \text{Dirichlet}(k \cdot \mathbf{w}_{\text{default}})$$
   Perturbs weights locally around default parameters ($\sim 95\%$ of draws within $\pm 0.10$ of baseline).
   - `python projects/NeuralTF/scripts/dirichlet_centered_all249.py` (all 289 candidates)
   - `python projects/NeuralTF/scripts/dirichlet_prioritize.py` (102 neural candidates)

2. **Uniform Dirichlet ($\alpha_i = 1$)**:
   $$\mathbf{w}^{(m)} \sim \text{Dirichlet}(\mathbf{1}_8)$$
   Samples uniformly across the entire 8-simplex to discover robust data-driven signals without prior preference.
   - `python projects/NeuralTF/scripts/dirichlet_uniform_all249.py` (all 289 candidates)
   - `python projects/NeuralTF/scripts/dirichlet_uniform.py` (102 neural candidates)

---

## ANANSE Gene Regulatory Network Scan

Validates candidates against the Perez 2025 ANANSE computational GRN (13,746 interactions across 9 cell fates):
- **30 / 289 candidates** act as primary upstream TF regulators
- **50 / 289 candidates** are downstream target genes in neural/differentiation pathways
- **12 candidates** function as core feedback hubs (both TF regulator and target)
- Output: `projects/NeuralTF/results/ananse_network_full.csv`

---

## Repository Structure

```
Bioinformatics/
├── pyproject.toml                            Package configuration & dependencies
├── README.md                                 Primary documentation
├── bioforge.md                               BioForge framework reference & operations
│
├── src/bioforge/                             BioForge Core Framework
│   ├── evidence/                             Multi-stream evidence engine
│   │   ├── schema.py                         EvidenceRecord & 8-stream EvidenceSource enum
│   │   ├── scoring.py                        Weighted score integration & DEFAULT_WEIGHTS
│   │   ├── confidence.py                     Tier classification (HIGH/MEDIUM/LOW)
│   │   └── cards.py                          Markdown evidence card generation
│   ├── projects/neuraltf/
│   │   ├── pipeline.py                       NeuralTFPipeline (5-atlas loader & 6 checkpoints)
│   │   ├── planmine.py                       PlanMine InterMine REST client & annotation parser
│   │   └── prioritize.py                     Dual-track candidate scoring & filtering
│   ├── omics/                                scRNA-seq QC, log-norm, PCA, Leiden clustering
│   ├── workflow/                             Declarative YAML workflow engine
│   ├── cli/                                  Command-line interface
│   └── ui/                                   Streamlit interactive dashboard
│
├── datasets/
│   ├── MANIFEST.md                           Download URLs & SHA256 checksums
│   ├── raw_data_manifest.md                  Complete 21-file raw data audit
│   ├── raw/                                  Raw downloads (gitignored)
│   └── processed/                            Processed H5ADs & Parquet files (gitignored)
│
├── projects/NeuralTF/
│   ├── data/
│   │   ├── bridge.csv                        v4 ↔ v6 ↔ gene_name Rosetta Stone mapping
│   │   ├── king_atlas.tsv                    King 2024 G0 progenitor enrichment table
│   │   ├── perez_tf_summary.csv              Perez 2025 TF lineage classification
│   │   └── master_tf_catalog.csv             Unified King + Perez master TF catalog (14,682 TFs)
│   ├── scripts/
│   │   ├── convert_cui.py                    Cui SMED → v6 H5AD (optional --subsample N)
│   │   ├── preprocess_perez.py               Parse Perez MOESM5 TF classes
│   │   ├── dirichlet_prioritize.py           Centered Dirichlet k=40 (102 neural TFs)
│   │   ├── dirichlet_uniform.py              Uniform Dirichlet α=1 (102 neural TFs)
│   │   ├── dirichlet_centered_all249.py      Centered Dirichlet k=40 (all 289 candidates)
│   │   ├── dirichlet_uniform_all249.py       Uniform Dirichlet α=1 (all 289 candidates)
│   │   ├── ananse_full_scan.py               ANANSE GRN scan across all 289 candidates
│   │   ├── export_fstf_ranked.py             Export ranked FSTF tables
│   │   ├── create_supplementary_tables.py    Generate supplementary tables S1–S5
│   │   ├── generate_publication_figures.py   Generate 21 publication figures
│   │   └── figures/                          21 modular figure generation scripts & style.py
│   ├── results/                              Dirichlet, ANANSE, and supplementary tables (gitignored)
│   ├── figures/                              21 publication-ready PNG figures (gitignored)
│   └── runs/pipeline_run/                    rank.csv, rank_neural.csv, 6 checkpoint parquets
│
└── scripts/                                  Master Orchestration & Build Scripts
    ├── generate_all.py                       End-to-end multi-step master pipeline runner
    ├── run_downstream.py                     Post-pipeline runner (Dirichlet, ANANSE, figures)
    ├── run.py                                Core pipeline execution entry point
    ├── build_bridge.py                       Build v4↔v6 gene ID bridge from Rosetta Stone
    ├── build_king_atlas.py                   Build king_atlas.tsv from King mmc7
    ├── build_master_catalog.py               Merge King mmc4 + Perez MOESM5 TF catalog
    ├── convert_fincher.py                    Convert Fincher DGE to H5AD
    └── consolidate_plass.py                  Consolidate Plass RAW.tar to H5AD
```

---

## Reproducibility & Verification

- **Deterministic Seeds**: All Dirichlet simulations run with `seed=2024` for 1,000 draws.
- **Fail-Fast Checkpoints**: Six Parquet checkpoints are recorded at each pipeline stage (`checkpoint_01` through `checkpoint_06`) in `projects/NeuralTF/runs/pipeline_run/`.
- **Unit Tests**: Full test suite passes:
  ```bash
  python -m pytest tests/
  ```

---

## License

MIT — see [LICENSE](LICENSE) file.
