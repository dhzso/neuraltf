# BioForge · NeuralTF

A reproducible pipeline for **planarian neural-fate-specific transcription factor** discovery. Integrates **five single-cell and regulatory atlases** (Fincher 2018, Plass 2018, Cui 2023, King 2024, Perez 2025) with 9 evidence streams, Bayesian Dirichlet uncertainty quantification, and ANANSE gene regulatory network validation to prioritize high-confidence targets for RNAi and functional validation.

> **2026-09 hardening audit**: the three prioritization methods (fixed /
> centered / uniform) now share one philosophy — the same all-candidate
> universe (rank.csv), the same annotation mask (one row per gene), the
> same +0.07 bonus layer (neural GO +0.03, TF GO +0.02, human ortholog
> +0.02), and the same Track-B DNA-binding-domain gate. The statistical
> suite reports circular AND circularity-controlled evaluations, and every
> figure is audited for non-emptiness
> (`projects/NeuralTF/scripts/audit_figures.py`).

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
python scripts/run_downstream.py             # Dirichlet UQ, ANANSE scan, tables & 33 figures

# 4. Launch the interactive Streamlit UI
bioforge ui                                  # http://localhost:8501
```

### Pipeline Outputs (`projects/NeuralTF/runs/pipeline_run/`)

| File | Content |
|------|---------|
| `rank.csv` | All **278 candidates** with scores across all 9 evidence streams |
| `rank_neural.csv` | **101 neural-enriched candidates** with proof status |
| `evidence_cards.md` | Per-candidate markdown evidence summary (278 cards) |
| `pipeline_results.json` | Machine-readable candidate metadata |
| `checkpoint_0*.parquet` | 6 incremental audit checkpoints (atlas loads, post-QC, post-scoring, King, Perez, stream matrix) |

### Prioritization & Analysis Outputs (`projects/NeuralTF/results/`)

| File | Content |
|------|---------|
| `master_tf_catalog.csv` | Unified King + Perez TF catalog (14,682 unique v6 TFs) — in `projects/NeuralTF/data/` |
| `dirichlet_centered_full_rank.csv` | Full Dirichlet-centered (k=40) composite rank (all candidates, 1 row/gene) |
| `dirichlet_uniform_full_rank.csv` | Full Dirichlet-uniform (α=1) composite rank (all candidates, 1 row/gene) |
| `dirichlet_centered_top10.csv` | Dual-track top-10 shortlist under centered Dirichlet (5 Track A + 5 Track B) |
| `dirichlet_uniform_top10.csv` | Dual-track top-10 shortlist under uniform Dirichlet (5 Track A + 5 Track B) |
| `dirichlet_*_draw_scores.csv` | Per-candidate draw-score matrices (bootstrap CIs, convergence) |
| `ananse_network_full.csv` | ANANSE GRN scan across all candidates & 9 cell fates (RBH-mapped, neuron-share normalized) |
| `ananse_top_regulators.csv` | Top planarian neural regulators (neuron-fate out-degree first) |
| `tf_ranked_neural_top19.csv` | Top 19 TFs: neural-filtered candidates |
| `tf_ranked_all_top43.csv` | Top 43 TFs: all expression-filtered candidates |
| `tf_ranked_catalog_top74.csv` | Top 74 TFs: full King mmc4 catalog |
| `supplementary_table_S1` – `S7` | Method comparison, and fixed/centered/uniform rank tables across all candidates |
| `de_pvalues.parquet` | Per-gene best DE p per atlas (pipeline checkpoint; drives meta-analysis) |

---

## What It Does

The pipeline seeds candidate TFs across five single-cell and regulatory atlases, computes a 9-stream multi-evidence matrix, evaluates scoring stability under Dirichlet uncertainty sampling, and maps candidates into cell-fate regulatory circuits.

### Atlases Integrated (5)

| Atlas | Year | Modality / Data Type | Role in Pipeline |
|-------|------|----------------------|------------------|
| **Fincher** | 2018 | scRNA-seq (50.5K cells, v4 IDs) | Whole-animal cell-type expression & neural specificity ([GSE111764](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111764)) |
| **Plass** | 2018 | scRNA-seq (37.5K cells, v6 IDs) | Independent whole-anatomy replication & X1 dynamics ([GSE103633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103633)) |
| **Cui** | 2023 | scRNA-seq (55.0K cells, 8 time points) | High-resolution regeneration time-course expression ([OMIX](https://ngdc.cncb.ac.cn/omix/release/OMIX003867)) |
| **King** | 2024 | Single-cell TF catalog + RNAi screen | G0/X1 cluster enrichment, RNAi phenotypes (mmc5), TF pair correlations (mmc6) ([Cell Reports](https://www.sciencedirect.com/science/article/pii/S2211124724001712)) |
| **Perez** | 2025 | Lineage atlas + ANANSE GRNs | Lineage TF classification (MOESM5), ANANSE regulatory influence (MOESM19), and ANANSE GRN (MOESM22) ([Nature Communications](https://www.nature.com/articles/s41467-025-65712-0#Sec94)) |

---

## Evidence Streams & Scoring Model (9 Streams)

Scoring utilizes a transparent, weighted multi-evidence integration model. Weights renormalize over streams present for each candidate:

$$\text{Integrated Score} = \sum_{i \in \text{Present}} w_i \cdot s_i \Bigg/ \sum_{i \in \text{Present}} w_i$$

| # | Stream | Default Weight ($w_i$) | Biological Basis & Computation |
|---|--------|------------------------|---------------------------------|
| 1 | **Expression** | 0.200 | $\min(1.0, \max(\text{log}_2\text{FC})/5)$ across Fincher, Plass, Cui, and King scRNA-seq atlases |
| 2 | **Specificity** | 0.100 | $1 / n_{\text{clusters}}$ supporting differential expression |
| 3 | **Reproducibility** | 0.100 | $n_{\text{atlases supporting}} / 5$ (Fincher, Plass, Cui, King, Perez) |
| 4 | **RNAi** | 0.100 | Binary indicator (1.0) if functional phenotype observed in King mmc5 screen |
| 5 | **Correlation** | 0.100 | $\min(1.0, \Delta r_{\text{G0-X1}} \times 3.0)$ co-expression gain from King mmc6 |
| 6 | **Neural Enriched** | 0.100 | Binary indicator (1.0) for G0 neural subcluster log₂FC ≥ 2.0 (King mmc7) |
| 7 | **Neural Specificity** | 0.100 | $1 / n_{\text{neural subclusters}}$ present in King atlas |
| 8 | **Perez Lineage** | 0.100 | Perez 2025 lineage TF class: **1.0** for neural-class, **0.5** for other TF classes, **0.0** if absent |
| 9 | **Perez Influence** | 0.100 | Perez 2025 ANANSE regulatory influence in neuron fate (MOESM19), normalized 0–1 rank |

### Confidence Tiers & Proof Status

- **Tier Assignment**:
  - **HIGH**: RNAi-validated OR (supporting streams ≥ 3 AND score ≥ 0.45)
  - **MEDIUM**: supporting streams ≥ 2 AND score ≥ 0.25
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
| 2 | `dd_Smed_v6_38342_0_1` | **dd38342** | 0.782 | 0.783 | `known_rnai_validated` | POU / Homeobox |
| 3 | `dd_Smed_v6_14115_0_1` | **dd14115** | 0.764 | 0.762 | `known_rnai_validated` | Homeobox / LIM |
| 4 | `dd_Smed_v6_14824_0_1` | **dd14824** | 0.756 | 0.762 | `known_rnai_validated` | C2H2 Zinc Finger |
| 5 | `dd_Smed_v6_19890_0_1` | **dd19890** | 0.753 | 0.753 | `known_rnai_validated` | T-box |
| 6 | `dd_Smed_v6_31217_0_1` | **dd31217** | 0.751 | 0.760 | `novel_candidate` | bHLH |
| 7 | `dd_Smed_v6_7033_0_1` | **dd7033** | 0.734 | 0.734 | `novel_candidate` | Homeobox / C2H2 ZNF |
| 8 | `dd_Smed_v6_4048_0_1` | **dd4048** | 0.750 | 0.749 | `novel_candidate` | bHLH |
| 9 | `dd_Smed_v6_11930_0_1` | **dd11930** | 0.734 | 0.734 | `novel_candidate` | C2H2 Zinc Finger |
| 10 | `dd_Smed_v6_9596_0_1` | **dd9596** | 0.691 | 0.691 | `novel_candidate` | Homeobox |

---

## Uncertainty Quantification & Sensitivity (Dirichlet Sampling)

To test sensitivity against arbitrary weighting assumptions, we employ Monte Carlo Dirichlet weight perturbation (1,000 draws, seed=2024 across all candidates in `rank.csv`):

1. **Centered Dirichlet ($k = 40$)**:
   $$\mathbf{w}^{(m)} \sim \text{Dirichlet}(k \cdot \mathbf{w}_{\text{default}})$$
   Perturbs weights locally around default parameters ($\sim 95\%$ of draws within $\pm 0.10$ of baseline).
   - `python projects/NeuralTF/scripts/dirichlet_centered.py`

2. **Uniform Dirichlet ($\alpha_i = 1$)**:
   $$\mathbf{w}^{(m)} \sim \text{Dirichlet}(\mathbf{1}_9)$$
   Samples uniformly across the entire 9-simplex to discover robust data-driven signals without prior preference.
   - `python projects/NeuralTF/scripts/dirichlet_uniform.py`

---

## Statistical Validation Suite (14 Tests)

The pipeline includes a comprehensive statistical validation suite to ensure publication-grade rigor:

| # | Test | Script | Key Metrics |
|---|------|--------|-------------|
| 1 | **Full Permutation Test** | `scripts/stats/permutation_test_full.py` | Empirical p-values (n=1,000 permutations) |
| 2 | **Bootstrap Confidence Intervals** | `scripts/stats/bootstrap_confidence.py` | 95% CI on integrated scores |
| 3 | **Overlap Significance** | `scripts/stats/overlap_significance.py` | Hypergeometric, Fisher's exact, binomial tests |
| 4 | **Precision-Recall Analysis** | `scripts/stats/precision_recall.py` | Precision@5, Precision@10, PR-AUC |
| 5 | **Negative Controls** | `scripts/stats/negative_controls.py` | Random non-TF & non-neural TF distributions |
| 6 | **Effect Sizes** | `scripts/stats/effect_sizes.py` | Cliff's delta, Cohen's d, Mann-Whitney U |
| 7 | **Leave-One-Atlas-Out** | `scripts/stats/leave_one_atlas_out.py` | Top-10 stability per excluded atlas |
| 8 | **Meta-Analytic P-values** | `scripts/stats/meta_analytic_pvalue.py` | Fisher's & Stouffer's combined p-values |
| 9 | **Power Analysis** | `scripts/stats/power_analysis.py` | Convergence & power curves |
| 10 | **Mann-Whitney U (Top-10)** | `scripts/stats/mann_whitney_top10.py` | Rank-biserial correlation |
| 11 | **Calibration Analysis** | `scripts/stats/calibration.py` | Empirical positive rates per decile |
| 12 | **Brier Score** | `scripts/stats/brier_score.py` | Probabilistic classification accuracy |
| 13 | **Cross-Method Correction** | `scripts/stats/cross_method_correction.py` | Bonferroni/BH-FDR for 3-method consensus |
| 14 | **Score Shuffling Permutation** | `scripts/stats/score_shuffling_permutation.py` | Stream-assignment null model |

Run all tests:
```bash
python scripts/run_statistical_tests.py
```

---

## Publication Figures (33 Figures)

### Main Figures (1–21)

| # | File | Description |
|---|------|-------------|
| 01 | `01_stream_coverage_all.png` | Evidence stream coverage across all TF candidates |
| 02 | `02_integrated_vs_composite.png` | Integrated vs composite score distribution |
| 03 | `03_score_distribution_all_vs_neural.png` | Score distribution: all vs neural-filtered |
| 04 | `04_evidence_heatmap_neural.png` | Evidence heatmap for neural candidates |
| 05 | `05_top10_candidate_atlas.png` | Top-10 candidate atlas visualization |
| 06 | `06_weight_sensitivity_ranks.png` | Weight sensitivity rank distributions |
| 07 | `07_weight_sensitivity_ptop10.png` | Weight sensitivity P(Top10) |
| 08 | `08_stream_ablation_global.png` | Stream ablation global impact |
| 09 | `09_stream_ablation_candidate.png` | Stream ablation candidate sensitivity |
| 10 | `10_centered_top10_scores.png` | Centered Dirichlet top-10 scores |
| 11 | `11_centered_scatter_neural.png` | Fixed vs centered Dirichlet scatter |
| 12 | `12_uniform_top10_scores.png` | Uniform Dirichlet top-10 scores |
| 13 | `13_uniform_scatter_all.png` | Fixed vs uniform Dirichlet scatter |
| 14 | `14_uniform_neural_vs_all_rankrank.png` | Neural vs all rank-rank comparison |
| 15 | `15_method_bumpchart.png` | 3-method rank comparison |
| 16 | `16_method_score_density.png` | 3-method score density |
| 17 | `17_method_rank_correlation.png` | 3-method rank correlation |
| 18 | `18_composite_bonus_waterfall.png` | Composite bonus waterfall |
| 19 | `19_method_consensus.png` | Method consensus |
| 20 | `20_stream_correlation.png` | Stream correlation matrix |
| 21 | `21_centered_vs_uniform_scatter.png` | Centered vs uniform Dirichlet |

### Statistical Figures (22–33)

| # | File | Description |
|---|------|-------------|
| 22 | `22_pipeline_schematic.png` | Conceptual pipeline diagram |
| 23 | `23_roc_pr_curve.png` | ROC and PR curves (RNAi ground truth) |
| 24 | `24_negative_controls.png` | Negative control distributions |
| 25 | `25_bootstrap_ci.png` | Bootstrap confidence intervals |
| 26 | `26_permutation_null.png` | Permutation null distribution |
| 27 | `27_loo_atlas_stability.png` | Leave-one-atlas-out stability |
| 28 | `28_effect_sizes.png` | Effect size annotations |
| 29 | `29_convergence_analysis.png` | Convergence analysis |
| 30 | `30_calibration.png` | Calibration reliability diagram |
| 31 | `31_score_distribution_all9.png` | Score distribution (9 streams) |
| 32 | `32_perez_influence_comparison.png` | Perez influence comparison |
| 33 | `33_method_agreement_summary.png` | Method agreement summary |

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
│   │   ├── schema.py                         EvidenceRecord & 9-stream EvidenceSource enum
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
│   │   ├── dirichlet_centered.py             Centered Dirichlet k=40 (all candidates)
│   │   ├── dirichlet_uniform.py              Uniform Dirichlet α=1 (all candidates)
│   │   ├── ananse_full_scan.py               ANANSE GRN scan across all candidates
│   │   ├── export_fstf_ranked.py             Export ranked TF tables
│   │   ├── create_supplementary_tables.py    Generate supplementary tables S1–S4
│   │   ├── generate_publication_figures.py   Generate 33 publication figures
│   │   └── figures/                          33 modular figure generation scripts & style.py
│   ├── results/                              Dirichlet, ANANSE, and supplementary tables (gitignored)
│   ├── figures/                              33 publication-ready PNG figures (gitignored)
│   └── runs/pipeline_run/                    rank.csv, rank_neural.csv, 6 checkpoint parquets
│
└── scripts/                                  Master Orchestration & Build Scripts
    ├── generate_all.py                       End-to-end multi-step master pipeline runner
    ├── run_downstream.py                     Post-pipeline runner (Dirichlet, ANANSE, figures, stats)
    ├── run_statistical_tests.py              Run all 14 statistical tests
    ├── run.py                                Core pipeline execution entry point
    ├── build_bridge.py                       Build v4↔v6 gene ID bridge from Rosetta Stone
    ├── build_king_atlas.py                   Build king_atlas.tsv from King mmc7
    ├── build_master_catalog.py               Merge King mmc4 + Perez MOESM5 TF catalog
    ├── convert_fincher.py                    Convert Fincher DGE to H5AD
    ├── consolidate_plass.py                  Consolidate Plass RAW.tar to H5AD
    └── stats/                                14 statistical test scripts
        ├── permutation_test_full.py
        ├── bootstrap_confidence.py
        ├── overlap_significance.py
        ├── precision_recall.py
        ├── negative_controls.py
        ├── effect_sizes.py
        ├── leave_one_atlas_out.py
        ├── meta_analytic_pvalue.py
        ├── power_analysis.py
        ├── mann_whitney_top10.py
        ├── calibration.py
        ├── brier_score.py
        ├── cross_method_correction.py
        └── score_shuffling_permutation.py
```

---

## License

MIT — see [LICENSE](LICENSE) file.
