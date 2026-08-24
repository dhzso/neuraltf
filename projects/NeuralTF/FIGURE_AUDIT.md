# Figure Audit — NeuralTF Visualization Redesign

## Methodology

Every existing figure was audited for:
1. **Scientific justification**: Does it answer a unique question?
2. **Statistical correctness**: Are scales, identity lines, and comparisons valid?
3. **Redundancy**: Does it duplicate information from another figure?
4. **Visual quality**: Is it publication-ready for Nature/Cell/Science family?

---

## Audit Results

### Old Figures from `visualize_fixed.py` (13 figures)

| # | Old Filename | Decision | Reason |
|---|---|---|---|
| 1 | `fig_fixed_score_distributions.png` | **REMOVED** | Histograms of 7 streams + integrated are low information density. Replaced by Fig1B (ECDF) and Fig2A (evidence heatmap). |
| 2 | `fig_fixed_candidate_summary.png` | **REMOVED** | 2×2 bar/box/line is a generic summary. Information is captured more precisely in Fig1 (landscape) and Fig2 (evidence architecture). |
| 3 | `fig_fixed_top10_dual_track.png` | **REMOVED** | Simple bar chart of Track A vs B scores. Captured in Fig6B (cross-method comparison) with more information. |
| 4 | `fig_fixed_evidence_heatmap.png` | **REPLACED** → Fig2A | Original showed all 99 candidates. New version focuses on Top 10 with metadata tracks and cell values. More informative at manuscript scale. |
| 5 | `fig_fixed_candidate_funnel.png` | **REPLACED** → Fig1A | Generic PowerPoint-style funnel replaced by publication-quality filtering flow with counts and percentages. |
| 6 | `fig_fixed_evidence_composition.png` | **REMOVED** | Stacked bars of evidence composition are redundant with Fig2B (weighted contribution decomposition) which is scientifically more meaningful. |
| 7 | `fig_fixed_stream_ablation.png` | **REPLACED** → Fig4 | Box+jitter replaced by lollipop plot (global impact) + candidate sensitivity heatmap. More statistically rigorous. |
| 8 | `fig_fixed_top10_radar.png` | **REMOVED** | Radar/polar charts are notoriously hard to read and compare. Information is better conveyed by Fig2A (heatmap) and Fig6B (lollipop). |
| 9 | `fig_fixed_go_dotplot.png` | **MOVED TO SUPPLEMENT** | GO term counts are secondary evidence. Moved to supplementary GO figures. |
| 10 | `fig_fixed_integrated_vs_composite.png` | **REMOVED** | Scatter of integrated vs composite is trivially expected (composite = integrated + small bonuses). Adds no scientific insight. |
| 11 | `fig_fixed_proof_status_violin.png` | **MOVED TO SUPPLEMENT** | Proof-status distributions are supplementary information, not a main figure question. |
| 12 | `fig_fixed_weight_sensitivity.png` | **REPLACED** → Fig3B | Line plot of rank trajectories replaced by rank uncertainty intervals with P(Top10) annotations. Statistically stronger. |
| 13 | `fig_fixed_integrated_vs_neural_filter.png` | **MOVED TO SUPPLEMENT** | ECDF comparison of 249 vs 99 is captured more comprehensively in Fig5B. |

### Old Figures from `visualize_centered.py` (5 figures)

| # | Old Filename | Decision | Reason |
|---|---|---|---|
| 1 | `fig_centered_trackA_top5.png` | **REMOVED** | Domain-colored bar chart is captured in Fig6B (cross-method comparison). |
| 2 | `fig_centered_trackB_top5.png` | **REMOVED** | Same as above. |
| 3 | `fig_centered_scatter_fixed_vs_dirichlet.png` | **REMOVED** | Scatter is captured in Fig3A (bump chart) which is more informative for rank comparison. |
| 4 | `fig_centered_combined_dual_track.png` | **REMOVED** | Vertical bar is redundant with Fig6B. |
| 5 | `fig_centered_score_shift.png` | **REMOVED** | Score shifts are not the right comparison; rank shifts are. Captured in Fig3. |

### Old Figures from `visualize_uniform.py` (7 figures)

| # | Old Filename | Decision | Reason |
|---|---|---|---|
| 1 | `fig_uniform_trackA_top5.png` | **REMOVED** | Redundant with Fig6B. |
| 2 | `fig_uniform_trackB_top5.png` | **REMOVED** | Redundant with Fig6B. |
| 3 | `fig_uniform_scatter_fixed_vs_uniform.png` | **REMOVED** | Captured in Fig3A and Fig5A. |
| 4 | `fig_uniform_scatter_centered_vs_uniform.png` | **REMOVED** | Captured in Fig3A. |
| 5 | `fig_uniform_combined_dual_track.png` | **REMOVED** | Redundant with Fig6B. |
| 6 | `fig_uniform_score_shift.png` | **REMOVED** | Rank shifts are more meaningful. Captured in Fig3. |
| 7 | `fig_uniform_three_way_comparison.png` | **REMOVED** | Grouped bar chart is replaced by Fig6B (lollipop) which is more readable. |

### Old Figures from `visualize_method_comparison.py` (5 figures)

| # | Old Filename | Decision | Reason |
|---|---|---|---|
| 1 | `fig_method_score_density.png` | **MOVED TO SUPPLEMENT** | KDE of score distributions is supplementary to the main ranking story. |
| 2 | `fig_method_rank_correlation.png` | **MOVED TO SUPPLEMENT** | Correlation matrices are supplementary statistical detail. |
| 3 | `fig_method_score_volatility.png` | **REMOVED** | Score volatility is not scientifically meaningful when methods use different weight distributions. |
| 4 | `fig_method_summary.png` | **REMOVED** | Overcrowded 2×3 subplot attempts to show everything. Each panel is now covered by a dedicated figure. |
| 5 | `fig_method_99vs249.png` | **REPLACED** → Fig5 | Simple bar+scatter replaced by rank-rank comparison + ECDF with statistical tests. |

### Old Supplementary Figures (4 figures)

| # | Old Filename | Decision | Reason |
|---|---|---|---|
| 1 | `fig_s1_go_gene_term_map.png` | **KEPT** | Well-designed binary heatmap. Scientifically valid. |
| 2 | `fig_s2_go_top10_dotmatrix.png` | **KEPT** | Focused GO matrix for top candidates. |
| 3 | `fig_s3_go_top_terms.png` | **KEPT** | Bar charts of top GO terms. |
| 4 | `fig_s4_go_neural_focus.png` | **KEPT** | Neural-focused GO heatmap. |

---

## Summary of Changes

| Category | Count |
|---|---|
| **Total old figures** | 34 |
| **Removed (redundant/low-value)** | 22 |
| **Replaced with new figures** | 6 |
| **Moved to supplement** | 4 |
| **Kept as-is** | 4 |
| **New publication figures** | 6 |
| **Net reduction** | 34 → 10 (6 main + 4 supplement) |

---

## Scoring Issues Discovered

1. **`rank.csv` and `rank_neural.csv` are identical** — both contain 99 rows. The pipeline does not produce a separate 249-candidate file. The full 249 candidates must be reconstructed from the pipeline's internal scoring.

2. **`correlation` stream is ALL NaN** — 6 of 7 evidence streams contribute. This means the correlation weight (0.105) is always redistributed to other streams via renormalization.

3. **`fixed_weight_score` not in `rank.csv`** — Only appears in Dirichlet comparison CSVs. The pipeline outputs `integrated_score` which is the same quantity.

4. **Dirichlet median scores are NOT on the same scale as fixed-weight scores** — While both are nominally 0–1, they use different weight vectors. Direct scatterplot comparison with y=x lines is valid only for composite scores (same formula applied).

5. **`dirichlet_uniform.py` line 226 bug** — Calls `read_mmc4()` instead of `read_mmc5()` for mmc5 loading. This may produce incorrect results for the mmc5 data integration step.

6. **`create_supplementary_tables.py` line 90 bug** — Tries to read a `.txt` file as CSV.
