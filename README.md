# BioForge · NeuralTF

A reproducible pipeline for **planarian neural-fate-specific transcription factor** discovery. Integrates **five single-cell and regulatory atlases** (Fincher 2018, Plass 2018, Cui 2023, King 2024, Perez 2025) — including the previously unused Fincher *brain* sub-atlas and the Cui *regeneration time-course* — into **11 evidence streams**, with Bayesian Dirichlet uncertainty quantification, a cross-species ortholog benchmark, and ANANSE gene regulatory network validation to prioritize high-confidence targets for RNAi and functional validation.

> **2026-09-11 ground-truth correction (screened ≠ validated)**: King 2024
> mmc5 is titled **"All Transcription Factors Inhibited"** — it lists every
> TF for which RNAi was *performed* and markers assayed, not only TFs with
> phenotypes. The table's own caption encodes phenotype status by font
> colour ("Transcription factors in red and marker genes in green have
> phenotypes"), but the distributed xlsx copy is monochrome. The historical
> `tested` label therefore means **screened**, not RNAi-validated. The
> FISH-confirmed phenotype set (**20 TFs**, curated from the paper's Fig
> 3J/4E, S4, S7, S8) is now tracked explicitly via the
> `phenotype_confirmed` column on every rank/shortlist table and via
> `bioforge.evidence.groundtruth.PHENOTYPE_CONFIRMED_V6`. All
> label-based statistics (ROC-AUC, calibration, effect sizes, negative
> controls) now report **both** labels: `screened` (legacy continuity)
> and `phenotype_confirmed` (the publishable "validated recovery"
> numbers). Under the stricter label the headline discrimination claims
> hold and strengthen: honest ROC-AUC 0.875, strict 0.833 (vs
> 0.846/0.807 on the screened label). The mmc7 green-font "previously
> published FSTF" annotation is also now captured (`king_atlas.tsv`:
> `previously_reported_fstf`, `sig_color`).

> **2026-09 hardening audit**: the three prioritization methods (fixed /
> centered / uniform) now share one philosophy — the same all-candidate
> universe (rank.csv), the same annotation mask (one row per gene), the
> same +0.07 bonus layer (neural GO +0.03, TF GO +0.02, human ortholog
> +0.02), and the same Track-B DNA-binding-domain gate. The statistical
> suite reports circular AND circularity-controlled evaluations, and every
> figure is audited for non-emptiness
> (`projects/NeuralTF/scripts/audit_figures.py`).
>
> **2026-09 evidence model v2**: candidate testing status is now reported as
> **Tested / Not tested / Known FSTF** (replacing the misleading "Track A/B" /
> "novel" framing; † = status read from the source paper's own tables), and the
> 9-stream model is extended to **11 streams** — two independent streams from
> previously unused raw data (`fincher_brain`, `cui_temporal`) — with weights
> rebalanced to sum to 1.0 (expression, `rnai`, and `correlation` were
> down-weighted so no single stream dominates). An independent cross-species
> ortholog benchmark (`ortholog_benchmark.py`) is added.

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
python scripts/run_downstream.py             # Dirichlet UQ, ANANSE scan, tables & 25 figures

# 4. Launch the interactive Streamlit UI
bioforge ui                                  # http://localhost:8501
```

### Pipeline Outputs (`projects/NeuralTF/runs/pipeline_run/`)

| File | Content |
|------|---------|
| `rank.csv` | All **11,695 candidates** with scores across all 11 evidence streams + `phenotype_confirmed` ground-truth flag |
| `rank_neural.csv` | Neural-enriched candidates with proof status + `ground_truth_status` |
| `evidence_cards.md` | Per-candidate markdown evidence summary |
| `pipeline_results.json` | Machine-readable candidate metadata |
| `de_pvalues.parquet` | Per-gene best DE p per atlas (meta-analysis input) |
| `checkpoint_0*.parquet` | 6 incremental audit checkpoints (atlas loads, post-QC, post-scoring, King, Perez, stream matrix) |

### Prioritization & Analysis Outputs (`projects/NeuralTF/results/`)

| File | Content |
|------|---------|
| `master_tf_catalog.csv` | Unified King + Perez TF catalog (14,682 unique v6 TFs) — in `projects/NeuralTF/data/` |
| `dirichlet_centered_full_rank.csv` | Full Dirichlet-centered (k=40) composite rank (all candidates, 1 row/gene) |
| `dirichlet_uniform_full_rank.csv` | Full Dirichlet-uniform (α=1) composite rank (all candidates, 1 row/gene) |
| `dirichlet_centered_top10.csv` | Dual-track top-10 shortlist under centered Dirichlet (5 tested + 5 not-tested) |
| `dirichlet_uniform_top10.csv` | Dual-track top-10 shortlist under uniform Dirichlet (5 tested + 5 not-tested) |
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

The pipeline seeds candidate TFs across five single-cell and regulatory atlases, computes an 11-stream multi-evidence matrix, evaluates scoring stability under Dirichlet uncertainty sampling, and maps candidates into cell-fate regulatory circuits.

### Atlases Integrated (5)

| Atlas | Year | Modality / Data Type | Role in Pipeline |
|-------|------|----------------------|------------------|
| **Fincher** | 2018 | scRNA-seq (50.5K cells, v4 IDs) | Whole-animal cell-type expression & neural specificity ([GSE111764](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111764)) |
| **Plass** | 2018 | scRNA-seq (37.5K cells, v6 IDs) | Independent whole-anatomy replication & X1 dynamics ([GSE103633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103633)) |
| **Cui** | 2023 | scRNA-seq (55.0K cells, 8 time points) | High-resolution regeneration time-course expression ([OMIX](https://ngdc.cncb.ac.cn/omix/release/OMIX003867)) |
| **King** | 2024 | Single-cell TF catalog + RNAi screen | G0/X1 cluster enrichment, RNAi phenotypes (mmc5), TF pair correlations (mmc6) ([Cell Reports](https://www.sciencedirect.com/science/article/pii/S2211124724001712)) |
| **Perez** | 2025 | Lineage atlas + ANANSE GRNs | Lineage TF classification (MOESM5), ANANSE regulatory influence (MOESM19), and ANANSE GRN (MOESM22) ([Nature Communications](https://www.nature.com/articles/s41467-025-65712-0#Sec94)) |

---

## Evidence Streams & Scoring Model (11 Streams)

Scoring utilizes a transparent, weighted multi-evidence integration model. Weights renormalize over streams present for each candidate:

$$\text{Integrated Score} = \sum_{i \in \text{Present}} w_i \cdot s_i \Bigg/ \sum_{i \in \text{Present}} w_i$$

| # | Stream | Default Weight ($w_i$) | Biological Basis & Computation |
|---|--------|------------------------|---------------------------------|
| 1 | **Expression** | 0.100 | $\min(1.0, \max(\text{log}_2\text{FC})/5)$ across Fincher, Plass, Cui, and King scRNA-seq atlases |
| 2 | **Specificity** | 0.100 | $1 / n_{\text{clusters}}$ supporting differential expression |
| 3 | **Reproducibility** | 0.100 | $n_{\text{atlases supporting}} / 5$ (Fincher, Plass, Cui, King, Perez) |
| 4 | **RNAi** | 0.050 | Binary indicator (1.0) if the TF was RNAi-screened with markers assayed in King mmc5 ("All Transcription Factors Inhibited") — screened, phenotype NOT implied |
| 5 | **Correlation** | 0.050 | $\min(1.0, \Delta r_{\text{G0-X1}} \times 3.0)$ co-expression gain from King mmc6 |
| 6 | **Neural Enriched** | 0.100 | Binary indicator (1.0) for G0 neural subcluster log₂FC ≥ 1.5 (King mmc7) |
| 7 | **Neural Specificity** | 0.100 | $1 / n_{\text{neural subclusters}}$ present in King atlas |
| 8 | **Perez Lineage** | 0.100 | Perez 2025 lineage TF class: **1.0** for neural-class, **0.5** for other TF classes, **0.0** if absent |
| 9 | **Perez Influence** | 0.100 | Perez 2025 ANANSE regulatory influence in neuron fate (MOESM19), normalized 0–1 rank |
| 10 | **Fincher Brain** | 0.100 | Independent Fincher 2018 BrainClustering sub-atlas enrichment ($\min(1.0, \text{log}_2\text{FC}/5)$ over brain Leiden clusters) |
| 11 | **Cui Temporal** | 0.100 | Cui 2023 regeneration time-course: neuronal temporal induction $\min(1.0, \log_2(\text{peak}/\text{baseline})/2)$ over the 8 timepoints |

### Confidence Tiers & Proof Status

- **Tier Assignment**:
  - **HIGH**: RNAi-screened OR (supporting streams ≥ 3 AND score ≥ 0.45)
  - **MEDIUM**: supporting streams ≥ 2 AND score ≥ 0.25
  - **LOW**: All other candidates
- **Testing Status** († = status read from the source paper's own tables):
  - `tested` — RNAi performed in the King et al. screen (†). **Screened,
    phenotype NOT implied** (see the ground-truth correction above).
  - `phenotype_confirmed` (column, alongside) — FISH-confirmed
    loss-of-cell-type phenotype in the paper's figures (20 TFs) †
  - `not_tested` — no RNAi record in the screen; *untested* rather than
    "novel" (many carry conserved orthologs)
  - `known_fstf` — documented fate-specifying TF (FSTF) from literature, no RNAi data (†)

---

## Top Prioritized Candidates (Production Run — Full Atlases, 11,695 Candidates)

Consensus across all three unified methods (fixed / centered Dirichlet k=40 / uniform Dirichlet α=1; shared universe, bonus mask, and gates) under the 11-stream model. 9 genes appear in the top-10 of **all three** methods — note the three methods share most streams by design, so this consensus measures weight-robustness, not method independence (see `overlap_significance.json` caveat).

| Testing status | Consensus candidates | Evidence |
|:---:|:---|:---|
| **Tested #1** | `dd_Smed_v6_38342_0_1` (**dd38342** / POU3F4-class) | RNAi-screened (†); consistent top-1 across methods |
| **Tested #2** | `dd_Smed_v6_34144_0_1` (**dd34144** / TCF7-LEF) | RNAi-screened (†) |
| **Tested #3** | `dd_Smed_v6_29211_0_1` (**dd29211** / PRRX2) | RNAi-screened (†); **FISH phenotype-confirmed** (dd8060+ neuron loss) |
| **Tested #4** | `dd_Smed_v6_22163_0_1` (**dd22163** / UNCX) | RNAi-screened (†); **FISH phenotype-confirmed** (sert+ neuron loss) |
| **Tested #5** | `dd_Smed_v6_12722_0_1` (**dd12722** / BHLHE23) | RNAi-screened (†) |
| **Not tested #1** | `dd_Smed_v6_5882_0_1` (**dd5882**) | Homeodomain; consensus top not-tested |
| **Not tested #2** | `dd_Smed_v6_14362_0_1` (**dd14362** / PAX5) | Paired-domain; consensus |
| **Not tested #3** | `dd_Smed_v6_12170_0_1` (**dd12170** / FOXJ2) | Forkhead; consensus |
| **Not tested #4** | `dd_Smed_v6_16466_0_1` (**dd16466** / FOXG1) | Forkhead; oral/neural forebrain ortholog; 2/3 methods |
| **Not tested #5** | `dd_Smed_v6_2442_0_1` (**dd2442** / HOXC6) | Homeobox; consensus |

Per-method shortlists: `dirichlet_centered_top10.csv`, `dirichlet_uniform_top10.csv`, `top10_neural_tfs_prioritized.csv` (all carry `phenotype_confirmed`).

Ground-truth recovery under the 11-stream model (both labels, circularity-controlled — label-encoding streams excluded): **screened** label (67 genes): honest ROC-AUC **0.846**, strict **0.807**, vs 0.976 circular; **phenotype_confirmed** label (19 in-universe genes): honest ROC-AUC **0.875**, strict **0.833**. An independent cross-species ortholog benchmark (neural-fate vs non-neural human orthologs) gives ROC-AUC **0.712** (provisional, coverage-limited).

### Known limitations (honest reading list)

1. **"Tested" = screened.** 48 of 67 screened TFs have no published phenotype; only 20 genes are FISH phenotype-confirmed (`bioforge.evidence.groundtruth`). All "validated" claims must use the `phenotype_confirmed` flag.
2. **Recall is imperfect.** `dd_Smed_v6_48508_0_1` has a confirmed neuron-loss phenotype but never entered the candidate universe (no seeded evidence) — the funnel misses it.
3. **Cross-method consensus is not independence.** The three methods share 9 of 11 streams, the bonus mask, and the King floors; the hypergeometric overlap p-values are an upper bound on significance.
4. **`fincher_brain` presence, not value, carries most of its signal.** Its unconditional per-stream AUC (~0.41 against the screened label) is a missingness artifact: among genes where the stream is present, the value-only AUC is ~0.66. The per-stream table in `precision_recall.json` now reports both variants.5. **King-derived expression still enters the honest score.** The honest/strict AUCs exclude label-encoding streams but keep King mmc7 expression floors; they are an upper bound on true label-free discrimination.
6. **Ortholog benchmark is provisional** (~101 classified genes from PlanMine BLAST descriptions; a full Compara/DIOPT table is needed for a definitive benchmark).

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

## Publication Figures (25 active figures)

The authoritative, up-to-date catalog is
`projects/NeuralTF/figures/FIGURE_CATALOG.md` (25 single-panel 500-DPI figures
under the Nature Communications palette). Figure numbering reflects the active
set: data-integration (01, 03, 04, 20, 31), candidate atlas (05),
robustness/ablation (06, 07, 08, 09, 13, 15), scoring decomposition (18),
benchmark recovery (23, 24, 26, 27, 28, 30), convergence (29), lineage (32),
method agreement (33), GRN (34), cross-atlas meta-analysis (35), and
regeneration temporal dynamics (36). Figures 02, 10–12, 14, 16, 17, 19, 21, 22,
25 were retired during the single-panel refactor.

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
│   │   ├── schema.py                        EvidenceRecord & 11-stream EvidenceSource enum
│   │   ├── scoring.py                       Weighted score integration & DEFAULT_WEIGHTS
│   │   ├── confidence.py                    Tier classification (HIGH/MEDIUM/LOW)
│   │   ├── groundtruth.py                   King-2024 ground truth: screened vs phenotype-confirmed labels
│   │   └── cards.py                         Markdown evidence card generation
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
│   │   ├── annotate_phenotype_groundtruth.py  Stamp phenotype_confirmed on all outputs (idempotent)
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
