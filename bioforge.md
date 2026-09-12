# bioforge.md — Operations, Architecture & Multi-Atlas Deep Dive

Extended technical and operational reference for `BioForge` and the `NeuralTF` planarian (*Schmidtea mediterranea*) neural transcription factor prioritization engine. Quick start: see [README.md](README.md).

---

## 1. Operational Workflow

### 1.1 First-Time Setup

```bash
git clone https://github.com/dhzso/neuraltf.git
cd neuraltf
python -m venv .venv
.venv\Scripts\activate        # Windows PowerShell / CMD
# source .venv/bin/activate   # Linux / macOS
pip install -e ".[bio,streamlit]"
```

### 1.2 Installation Verification

```bash
# Smoke test core imports
python -c "from bioforge import *; print('BioForge OK')"

# Run core unit test suite
python -m pytest tests/unit/ -v -q
```

### 1.3 Preprocessing Raw Atlases

Raw dataset files are downloaded into `datasets/raw/` (see `datasets/MANIFEST.md` for download URLs and checksums). Build the processed H5AD and reference tables locally:

```bash
# 1. Fincher 2018 Atlas (GEO GSE111764) -> datasets/processed/fincher_subsample.h5ad
python scripts/convert_fincher.py

# 2. Plass 2018 Atlas (GEO GSE103633) -> datasets/processed/plass_v6.h5ad
python scripts/consolidate_plass.py

# 3. Cui 2023 Atlas (OMIX003867, full 55K cells) -> datasets/processed/cui_v6.h5ad
python projects/NeuralTF/scripts/convert_cui.py

# 4. Perez 2025 TF Lineage Classification (MOESM5) -> projects/NeuralTF/data/perez_tf_summary.csv
python projects/NeuralTF/scripts/preprocess_perez.py

# 5. Master TF Catalog Builder (King mmc4 + Perez MOESM5) -> projects/NeuralTF/data/master_tf_catalog.csv
python scripts/build_master_catalog.py
```

### 1.4 Running the Multi-Atlas NeuralTF Pipeline

Execute the unified 5-atlas prioritization pipeline:

```bash
bioforge neuraltf run [--subsample 0] [--out projects/NeuralTF/runs/pipeline_run]
# Or directly:
python scripts/run.py
```

#### Pipeline Outputs in `projects/NeuralTF/runs/pipeline_run/`

| File | Format | Description |
|------|--------|-------------|
| `rank.csv` | CSV (278 rows) | All **11,695 candidates** ranked by 11-stream integrated score |
| `rank_neural.csv` | CSV (101 rows) | **134 neural candidates** passing the neural gate |
| `evidence_cards.md` | Markdown | Comprehensive per-candidate evidence cards with stream breakdown |
| `pipeline_results.json` | JSON | Machine-readable candidate records with tier classifications |
| `checkpoint_01_atlas_loads.parquet` | Parquet | QC checkpoint: Atlas cell and gene dimensions |
| `checkpoint_02_post_qc.parquet` | Parquet | QC checkpoint: Post-filter and Leiden cluster assignments |
| `checkpoint_03_post_scoring.parquet` | Parquet | QC checkpoint: Per-atlas Wilcoxon DE scores |
| `checkpoint_04_king_records.parquet` | Parquet | QC checkpoint: King G0 progenitor neural subcluster records |
| `checkpoint_05_perez_records.parquet` | Parquet | QC checkpoint: Perez TF superfamily lineage scores |
| `checkpoint_06_stream_matrix.parquet` | Parquet | Full 11-stream evidence feature matrix across all 278 candidates |

---

### 1.5 Automated Downstream Analysis Pipeline

Run all downstream uncertainty quantification, network scans, supplementary tables, statistical tests, and publication figures in a single dependency-managed command:

```bash
python scripts/run_downstream.py
```

Or execute individual downstream modules:

```bash
# Centered Dirichlet uncertainty quantification (k=40 across all candidates)
python projects/NeuralTF/scripts/dirichlet_centered.py

# Uniform Dirichlet prior robustness scan (alpha=1 across all candidates)
python projects/NeuralTF/scripts/dirichlet_uniform.py

# ANANSE Gene Regulatory Network validation (Perez 2025 MOESM22)
python projects/NeuralTF/scripts/ananse_full_scan.py

# PlanMine functional annotation & dual-track prioritization
python scripts/prioritize_neural_tfs.py

# Statistical validation suite (14 tests)
python scripts/run_statistical_tests.py

# Supplementary tables & 33 publication figures
python projects/NeuralTF/scripts/create_supplementary_tables.py
python projects/NeuralTF/scripts/generate_publication_figures.py
```

---

### 1.6 Interactive Streamlit Dashboard

```bash
bioforge ui [--port 8501] [--host localhost]
```

Accessible at `http://localhost:8501` with four specialized tabs:
1. **Run Page**: Real-time dataset discovery, QC monitoring, and pipeline execution controls.
2. **Results Page**: Interactive rank tables, dynamic filtering, scatter/density plots, and evidence cards.
3. **Prioritization Page**: Dual-track candidate evaluation (Tested vs Not-tested candidates).
4. **AI Assistant**: Conversational biology assistant for hypothesis generation and candidate interpretation (Need *API_Key* config).

---

## 2. Multi-Atlas Architecture & Integration

NeuralTF unifies 5 independent planarian transcriptomic and regulatory atlases:

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                               5-ATLAS INPUT UNIVERSE                                   │
├──────────────────┬──────────────────┬──────────────────┬───────────────┬───────────────┤
│   Fincher 2018   │    Plass 2018    │     Cui 2023     │   King 2024   │  Perez 2025   │
│ (Science, Drop)  │ (Science, Drop)  │ (NatCom, 10x 55K)│ (CellRep, G0) │(NatCom, Linea)│
└────────┬─────────┴────────┬─────────┴────────┬─────────┴───────┬───────┴───────┬───────┘
         │                  │                  │                 │               │
         ▼                  ▼                  ▼                 ▼               ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                              11 EVIDENCE STREAMS MATRIX                                 │
├────────────────────────────────────────────────────────────────────────────────────────┤
│ 1. Expression ($w_1=0.100$)         : Best log2FC / 5.0 across scRNA-seq atlases       │
│ 2. Specificity ($w_2=0.100$)        : Inverse cluster breadth (1 / n_clusters)         │
│ 3. Reproducibility ($w_3=0.100$)    : Cross-atlas concordance (n_supporting / 5)       │
│ 4. RNAi ($w_4=0.050$)               : RNAi-screened in King mmc5 ("All TFs Inhibited")│
│ 5. Correlation ($w_5=0.050$)        : G0 vs X1 co-expression correlation gain         │
│ 6. Neural Enriched ($w_6=0.100$)    : King G0 neural subcluster log2FC ≥ 1.5 (1 / 0)   │
│ 7. Neural Specificity ($w_7=0.100$) : Inverse neural subcluster breadth (1 / n_subs)   │
│ 8. Perez Lineage ($w_8=0.100$)      : Perez TF structural class (1.0 / 0.5 / 0.0)      │
│ 9. Perez Influence ($w_9=0.100$)    : Perez ANANSE neuron-fate influence (0.0 to 1.0)  │
│ 10. Fincher Brain ($w_{10}=0.100$)  : Fincher 2018 BrainClustering enrichment (l2FC/5)│
│ 11. Cui Temporal ($w_{11}=0.100$)   : Cui 2023 neuronal temporal induction (peak/base)  │
└──────────────────────────────────────────┬─────────────────────────────────────────────┘
                                           │
                                           ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        COMPOSITE SCORING & UNCERTAINTY QUANTIFICATION                  │
├────────────────────────────────────────────────────────────────────────────────────────┤
│ • Baseline Integrated Score : Weighted linear combination (renormalized over present)  │
│ • Dirichlet UQ (Centered)   : 1,000 draws from Dirichlet(k=40 * w_default)             │
│ • Dirichlet UQ (Uniform)    : 1,000 draws from Dirichlet(alpha=1_9)                    │
│ • ANANSE GRN Validation     : 13,746 TF-target edges across 9 cell fate lineages       │
│ • Statistical Suite         : 14 validation tests (permutations, bootstrap, calibration)│
└──────────────────────────────────────────┬─────────────────────────────────────────────┘
                                           │
                                           ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                          DUAL-TRACK CANDIDATE PRIORITIZATION                           │
├──────────────────────────────────────────┬─────────────────────────────────────────────┤
│   Track A: Screened Benchmarks (Top 5)    │       Track B: Untested TFs (Top 5)       │
│   (RNAi-screened in King 2024 mmc5;      │       (High-scoring, domain-gated        │
│    phenotype-confirmed subset flagged †) │        discoveries)                       │
└──────────────────────────────────────────┴─────────────────────────────────────────────┘
```

---

## 3. Candidate Selection Funnel & Mathematical Formulation

### 3.1 Candidate Count Funnel

1. **418 TF Targets** (`load_reference_tables`): Seeded from King 2024 `mmc4.xlsx` TF catalog (`TF? != NA`).
2. **Cluster DE Candidates** (`score_atlases`): TFs displaying statistically significant differential expression (Wilcoxon rank-sum test, Benjamini-Hochberg $q \le 0.10$) across Leiden clusters in Fincher, Plass, and Cui's scRNA-seq atlases.
3. **11,695 Total Scored Candidates** (`integrate_king_atlas`): Integrating King 2024 `mmc7.xlsx` G0 neural subclusters ($\text{log}_2\text{FC} \ge 1.5$, matching King's own STAR Methods), King mmc5 RNAi screening targets, and Perez 2025 regulatory data seeds these factors into `all_records`.
4. **134 Neural Candidates** (`write_outputs`): Applying the neural filter `(neural_enriched > 0) | (rnai > 0)` yields **134 candidates** across Track A (RNAi-screened) and Track B (untested candidates).

---

### 3.2 Evidence Stream Mathematical Formulations

$$\text{Integrated Score}(g) = \frac{\sum_{i=1}^{11} w_i \cdot s_i(g) \cdot \mathbb{I}(s_i(g) \text{ present})}{\sum_{i=1}^{11} w_i \cdot \mathbb{I}(s_i(g) \text{ present})}$$

Where default weights $\mathbf{w} = [0.100, 0.100, 0.100, 0.050, 0.050, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100]$:

1. **Expression ($s_{\text{expr}}$)**:
   $$s_{\text{expr}} = \min\left(1.0, \frac{\max(\text{log}_2\text{FC})}{5.0}\right)$$

2. **Specificity ($s_{\text{spec}}$)**:
   $$s_{\text{spec}} = \frac{1.0}{n_{\text{clusters expressing}}}$$

3. **Reproducibility ($s_{\text{repro}}$)**:
   $$s_{\text{repro}} = \frac{\min(n_{\text{supporting atlases}}, 5)}{5.0}$$

4. **RNAi Phenotype ($s_{\text{rnai}}$)**:
   $$s_{\text{rnai}} = \mathbb{I}(g \in \text{King 2024 mmc5 RNAi screening table})$$
   *Ground-truth note (2026-09-11):* mmc5 is titled "All Transcription Factors Inhibited" — presence means RNAi was performed and markers assayed, NOT that a phenotype was observed (the original table encodes phenotype status by font colour; the distributed copy is monochrome). The FISH-confirmed subset is tracked via `bioforge.evidence.groundtruth.PHENOTYPE_CONFIRMED_V6` and the `phenotype_confirmed` output column.

5. **Co-Expression Correlation Gain ($s_{\text{corr}}$)**:
   $$s_{\text{corr}} = \min\left(1.0, \max(0.0, r_{\text{G0}} - r_{\text{X1}}) \times 3.0\right)$$

6. **Neural Enrichment ($s_{\text{neural-enr}}$)**:
   $$s_{\text{neural-enr}} = \mathbb{I}(\text{King G0 neural subcluster } \text{log}_2\text{FC} \ge 2.0)$$

7. **Neural Specificity ($s_{\text{neural-spec}}$)**:
   $$s_{\text{neural-spec}} = \frac{1.0}{n_{\text{neural subclusters with } \text{log}_2\text{FC} \ge 2.0}}$$

8. **Perez Lineage ($s_{\text{perez-lin}}$)**:
   - **1.0**: Assigned if TF class belongs to Neural Superfamilies (bHLH, Homeobox, POU, C2H2, SOX, FOX, etc.)
   - **0.5**: Assigned if TF class is another confirmed structural TF family
   - **0.0**: Assigned if unclassified or absent from Perez 2025 catalog

9. **Perez Regulatory Influence ($s_{\text{perez-infl}}$)**:
   - Evaluated as the normalized ANANSE regulatory influence score in neuron fate specification (range $[0.0, 1.0]$).
   - *Data source: [Perez et al. (2025)](https://www.nature.com/articles/s41467-025-65712-0#Sec94) MOESM19.*

---

## 4. Top Prioritized Candidates

Dual-track shortlist under the 11-stream model (updated 2026-09-11; see README for the canonical table). **Ground-truth note:** King mmc5 lists ALL inhibited TFs ("All Transcription Factors Inhibited"); the distributed xlsx lost its red/green phenotype font encoding, so `tested` means *screened*. The `phenotype_confirmed` column marks the 20 FISH-confirmed TFs (paper Fig 3J/4E, S4, S7, S8); only two of the five Track-A members below are phenotype-confirmed (dd29211/PRRX2, dd22163/UNCX).

| Rank | Gene ID | Gene Name | TF Class / Family | Integrated Score | Dirichlet Median ($k=40$) | Status | Phenotype-confirmed |
|:---:|:---|:---|:---|:---:|:---:|:---|:---:|
| **1** | `dd38342` | *dd38342* | POU / Homeobox | **0.905** | ~0.90 | `tested` | — |
| **2** | `dd34144` | *dd34144* | HMG / TCF-LEF | **0.902** | ~0.90 | `tested` | — |
| **3** | `dd29211` | *dd29211* | Homeobox / PRRX2 | **0.881** | ~0.88 | `tested` | † |
| **4** | `dd22163` | *dd22163* | Homeobox / UNCX | **0.865** | ~0.86 | `tested` | † |
| **5** | `dd12722` | *dd12722* | bHLH / BHLHE23 | **0.863** | ~0.86 | `tested` | — |
| **6** | `dd5882` | *dd5882* | Homeodomain | **0.892** | ~0.89 | `not_tested` | — |
| **7** | `dd14362` | *dd14362* | Paired / PAX5 | **0.819** | ~0.82 | `not_tested` | — |
| **8** | `dd12170` | *dd12170* | Forkhead / FOXJ2 | **0.841** | ~0.84 | `not_tested` | — |
| **9** | `dd16466` | *dd16466* | Forkhead / FOXG1 | **0.857** | ~0.86 | `not_tested` | — |
| **10** | `dd2442` | *dd2442* | Homeobox / HOXC6 | **0.834** | ~0.83 | `not_tested` | — |

*Integrated scores shown are from the fixed-weight method; Dirichlet medians are approximate ($k=40$).*

---

## 5. Software Architecture & Directory Layout

```
src/bioforge/
├── core/                  # Configuration, logging, exception hierarchies
├── evidence/              # 11-stream scoring engine, EvidenceScorer, EvidenceRecord
├── projects/neuraltf/     # Multi-atlas pipeline, PlanMine client, prioritization engine
├── omics/                 # Single-cell QC, normalization, clustering, Leiden algorithms
├── smapping/              # Cross-assembly identifier mapping (SMED ↔ v4 ↔ v6 ↔ h1SMcG)
├── ui/                    # Streamlit multi-page analytical dashboard
└── ai/                    # LLM assistant integrations (OpenAI / Ollama / Stub)

projects/NeuralTF/
├── data/                  # bridge.csv, king_atlas.tsv, master_tf_catalog.csv, perez_tf_summary.csv
├── results/               # Dirichlet CSVs, ANANSE network, top-10 prioritization, tables S1–S4
├── figures/               # 25 active single-panel 500 DPI figures
└── runs/pipeline_run/     # rank.csv, rank_neural.csv, evidence_cards.md, audit checkpoints 01–06
```
