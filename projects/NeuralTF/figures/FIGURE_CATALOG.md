# NeuralTF Publication Figure Catalog

This catalog provides an authoritative reference for all **25 publication figures** generated for the NeuralTF project. All figures strictly comply with *Nature Communications* formatting guidelines:
- **Resolution**: 500 DPI rasterized PNG.
- **Format**: Exactly one standalone figure per PNG file (no collages, no embedded composite sub-panels).
- **Typography**: Clean Arial / Helvetica sans-serif; panel titles 8.0 pt regular/medium, axis labels 7.0 pt regular, tick labels & annotations 6.0–6.2 pt.
- **Color Palette**: Nature-standard muted palette (Deep navy `#1B365D`, terracotta `#C25943`, warm amber `#D9822B`, soft gray `#78909C`).
- **Data Integrity**: Honest, non-circular statistical baselines and exact gene nomenclature (`dd_Smed_v6` cross-referenced to canonical gene symbols).

---

## Active Figure Master Table

| Figure File | Panel Title | Category | Statistical / Biological Role |
| :--- | :--- | :--- | :--- |
| `01_stream_coverage_all.png` | Empirical Completeness Across Nine Evidence Streams | Data Integration | Whole-transcriptome completeness ($N=11,675$) across all 9 evidence streams. |
| `03_score_distribution_all_vs_neural.png` | Integrated Evidence Score Separation: All TFs vs Neural Regulators | Score Validation | Empirical KDE separating background transcriptome vs prioritized neural TFs (K-S test $D=0.901, P < 10^{-30}$) with subtitle test reporting, vertical dashed medians, and generous headroom. |
| `04_evidence_heatmap_neural.png` | Multi-Stream Evidence Landscape Across All 134 Prioritized Neural Transcription Factors | Candidate Landscape | Stream-by-stream normalized evidence heatmap for all 134 neural TFs across all canonical 9 evidence streams, with neutral gray missing-value cells, track groupings, isoform disambiguation, and integrated score bars. |
| `05_top10_candidate_atlas.png` | Candidate Atlas: High-Confidence & Novel Neural Transcription Factor Regulators | Candidate Profiling | Dual-panel atlas: (a) 9-stream canonical evidence matrix with numerical cell values; (b) Two-tone stacked prioritization bars (base evidence in track color + composite bonus in amber gold), TF family badges, and human ortholog annotations. |
| `06_weight_sensitivity_ranks.png` | Prior Sensitivity: Top Candidate Rank Uncertainty | Robustness | Rank distributions across 1,000 Centered Dirichlet weight draws ($k=40$) for prioritized candidates and top challengers, with explicit selection threshold and stable anchors. |
| `07_weight_sensitivity_ptop10.png` | Centered Dirichlet Prior Sensitivity: Top-10 Inclusion Probability | Robustness | Candidate inclusion probability $P(\text{Top } 10)$ across 1,000 Centered Dirichlet draws ($k=40$) for 58 candidates and challengers, with majority ($\ge 50\%$) and high-confidence ($\ge 80\%$) consensus thresholds. |
| `08_stream_ablation_global.png` | Global Stream Ablation: Spearman Rank Degradation | Stream Importance | Impact on transcriptome-wide rankings when each evidence stream is omitted. |
| `09_stream_ablation_candidate.png` | Stream Dependency Matrix: Top-10 Score Deltas | Stream Importance | Individual candidate score drops ($\Delta \text{score}$) upon single-stream ablation. |
| `13_uniform_scatter_all.png` | Centered ($k=40$) vs Uniform ($\alpha=1$) Dirichlet Scores | Prior Invariance | Genome-wide concordance demonstrating invariance to uninformative flat priors ($r_s = 0.999$). |
| `15_method_bumpchart.png` | Candidate Rank Trajectories Across Dirichlet Weighting | Method Concordance | Dual-track bump chart tracking within-track rank trajectories across Fixed, Centered Dirichlet ($k=40$), and Uniform Dirichlet ($\alpha=1$) weight formulations for high-confidence and novel candidates. |
| `18_composite_bonus_waterfall.png` | Additive Bonus Contributions to Candidate Prioritization | Scoring Decomposition | Exact breakdown: Base score + GO Neural bonus + GO TF bonus + Human Ortholog bonus, with distinct vertical clearance between title and legend. |
| `20_stream_correlation.png` | Pairwise Correlation Matrix of Evidence Streams | Stream Independence | Pairwise Spearman correlation matrix confirming low-to-moderate collinearity between streams. |
| `23_roc_curve.png` | Receiver Operating Characteristic: Neural TF Recovery | Benchmark Recovery | ROC curves comparing circular benchmark, honest label-free, and strict label-free models. |
| `24_negative_controls.png` | Biological Specificity: Score Distribution vs Empirical Negative Controls | Specificity | Label-free score distributions (excluding circular streams) in prioritized neural candidates ($n=67$) vs matched negative control cohorts ($n=100$) with jitter points, medians, two-sided Mann–Whitney $U$ tests ($U_1 = 5,620, U_2 = 5,468, P < 10^{-11}$), and Cohen's $d$ effect sizes ($d = 1.34, 1.30$). |
| `26_permutation_null.png` | Permutation Null vs Observed Candidate Scores | Hypothesis Testing | Empirical null distribution under joint row permutations vs shaded observed candidate range (0.94–1.00) with subtitle permutation test reporting ($P \le 0.0323$). |
| `27_pr_curve.png` | Precision–Recall: Neural TF Recovery | Benchmark Recovery | Precision-recall curves with empirical background prevalence baseline and circularity control. |
| `28_effect_sizes.png` | Effect Sizes and Significance Across Candidate Cohorts | Effect Magnitude | Non-parametric Cliff's $\delta$ and parametric Hedges' $g$ quantifying separation across cohorts with exact two-sided Mann–Whitney test p-values ($P \le 2.2 \times 10^{-8}$). |
| `29_convergence_analysis.png` | Rank Stability Convergence Under Dirichlet Resampling | Convergence | Spearman rank correlation vs number of draws, proving convergence ($r_s \geq 0.99$) at $n=200$. |
| `30_calibration.png` | Rank Discrimination and Score Calibration Across Transcriptome Deciles | Score Calibration | RNAi validation rate across score deciles (D1–D10) with Wilson 95% CIs, genome-wide prevalence reference ($0.57\%$), and top-decile 9.5x enrichment callout with exact binomial test significance ($P = 2.0 \times 10^{-40}$). |
| `31_score_distribution_all9.png` | Empirical Distributions Across Nine Evidence Streams | Data Architecture | Horizontal violin plots showing distribution, median, and IQR of non-zero scores per stream. |
| `32_perez_influence_comparison.png` | Evidence Score Stratification Across Lineages | Lineage Specificity | Evidence score distributions across single-cell lineage classes from Perez et al. (Neural, Other Lineages, and Unclassified background with $n=300$ subsampled jitter points, Mann–Whitney $P < 10^{-15}$). |
| `33_method_agreement_summary.png` | Prioritization Method Agreement (Top 10 Candidates) | Method Concordance | Pairwise Jaccard similarity matrix across Fixed, Centered Dirichlet, and Uniform Dirichlet models with shared candidate counts and exact hypergeometric overlap test p-values ($P < 10^{-13}$) and 3-way consensus binomial test ($P = 1.8 \times 10^{-28}$). |
| `34_ananse_regulatory_network.png` | ANANSE Regulatory Network Target Specificity | GRN Topology | Scatter plot of total predicted gene targets vs neuron-specific targets from ANANSE binding networks. |
| `35_meta_analysis_concordance.png` | Cross-Atlas Neural Effect Size Concordance | Reproducibility | Concordance of neural $\log_2$ fold change between Fincher et al. and Plass et al. scRNA-seq atlases. |
| `36_regeneration_temporal_dynamics.png` | Regeneration Temporal Dynamics (Cui et al. 2023) | Longitudinal Biology | Standardized expression ($z$-score) trajectories across 8 regeneration stages (0d to 7d post-amputation). |

---

## Directory Organization Note

- All current publication figure PNG files are stored directly in `projects/NeuralTF/figures/`.
- No raw CSV data files, incomplete sketches, or multi-panel publication collages reside in this folder.
- All generating scripts are maintained in `projects/NeuralTF/scripts/figures/` and can be re-run in batch via `projects/NeuralTF/scripts/generate_publication_figures.py`.
