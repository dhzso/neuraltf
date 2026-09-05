# NeuralTF Pipeline — Parameter Scientific Rationale & Sensitivity Documentation

This document records the biological, statistical, and mathematical justifications for all parameters, thresholds, and scoring formulations used in the `NeuralTF` pipeline.

---

## 1. Multi-Atlas Evidence Stream Scoring Formulations

### 1.1 Expression Score Cap ($\text{Divisor} = 5.0$)

$$\text{Score}_{\text{expression}} = \min\left(1.0, \frac{\max(\text{log}_2\text{FC})}{5.0}\right)$$

- **Biological Rationale**: In planarian (*Schmidtea mediterranea*) single-cell RNA (scRNA) sequencing, a $\text{log}_2\text{FC} = 5.0$ corresponds to a $2^5 = 32$-fold of linear upregulation which is relative to background non-target cells. Sometimes, transcription factor (TF) steady-state transcript abundance in whole-animal or tissue clusters reaches saturation in target promoter occupancy well before 32-fold upregulation. Setting the divisor at 5.0 maps the physiological range $[0, 5.0]$ linearly to $[0, 1.0]$. This avoids outlier compression where a single extreme fold-change (e.g., $\text{log}_2\text{FC} = 12$) would compress biologically meaningful, moderately expressed lineage master regulators ($\text{log}_2\text{FC} \in [2.0, 4.0]$) into near-zero scores.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:413-416` (Fincher/Plass/Cui), `pipeline.py:482-484` (King)
- **Sensitivity Analysis**:
  | Divisor | Total Candidates | Top-10 Jaccard vs 5.0 | Biological Behavior |
  |:---:|:---:|:---:|:---|
  | 3.0 | 289 | 0.88 | Premature saturation at 8-fold FC; compresses high-expression dynamic range |
  | **5.0** | **289** | **1.00** | **Optimal linear dynamic range (Default)** |
  | 8.0 | 289 | 0.95 | Down-weights moderately expressed but functional TFs |
  | 10.0 | 289 | 0.90 | Biases towards high-copy transcripts at the expense of regulatory factors |


### 1.2 Multi-Atlas Reproducibility Denominator ($N_{\text{atlases}} = 5$)

$$\text{Score}_{\text{reproducibility}} = \frac{\min(n_{\text{supporting atlases}}, 5)}{5.0}$$

- **Biological Rationale**: True lineage specification drivers exhibit consistent transcriptional induction across independent experimental protocols, cell dissociation methods, and sequencing technologies. The pipeline interrogates all 5 single-cell/regulatory atlases:
  1. **Fincher 2018** (Drop-seq, whole-animal, `dd_Smed_v4` transcript model)
  2. **Plass 2018** (Drop-seq, whole-animal, `dd_Smed_v6` transcript model)
  3. **Cui 2023** (10x Genomics Chromium, regeneration single-cell atlas, 55,014 cells)
  4. **King 2024** (FACS-purified G0/X1 progenitor single-cell atlas & TF catalog)
  5. **Perez 2025** (Lineage Single-Cell Differentiation Atlas & TF Domain Family Classification)
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:663-672`
- **Scientific Impact**: Candidate TFs supported across all 5 atlases achieve $s_{\text{repro}} = 1.0$, while single-atlas detections receive $s_{\text{repro}} = 0.20$.


### 1.3 King Atlas Neural Subcluster Enrichment Threshold ($\text{log}_2\text{FC} \ge 2.0$)

$$\text{Neural Gate} = \mathbb{I}\left(\text{subcluster} \in \text{Neural} \land \text{log}_2\text{FC} \ge 2.0\right)$$

- **Biological Rationale**: In planarian adult stem cell differentiation (King et al., *Cell Reports* 2024), pluripotent neoblasts (X1) exit the cell cycle into post-mitotic committed progenitors (G0). Fate-specifying transcription factors (FSTFs) undergo sharp, switch-like transcriptional activation during this transition. A $\text{log}_2\text{FC} \ge 2.0$ represents a $\ge 4$-fold enrichment in neural G0 progenitors compared to general G0 populations.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:441`
- **Sensitivity Analysis**:
  | $\text{log}_2\text{FC}$ Cutoff | Neural Candidates | Top-10 Jaccard vs 2.0 | Statistical Assessment |
  |:---:|:---:|:---:|:---|
  | 1.0 | 148 | 0.70 | Permissive; includes broad lineage progenitors and false positives |
  | 1.5 | 122 | 0.82 | Intermediate; captures minor sublineage markers |
  | **2.0** | **102** | **1.00** | **Optimal balance of sensitivity and specificity (Default)** |
  | 2.5 | 78 | 0.76 | Overly stringent; excludes validated neural TFs with modest basal expression |
  | 3.0 | 60 | 0.62 | Severe loss of known RNAi-validated regulators |


### 1.4 Co-Expression Correlation Gain Multiplier ($\text{Multiplier} = 3.0$)

$$\text{Score}_{\text{correlation}} = \min\left(1.0, \max(0.0, r_{\text{G0}} - r_{\text{X1}}) \times 3.0\right)$$

- **Biological Rationale**: In pluripotent X1 neoblasts, neural TF pairs display negligible co-expression correlation ($r_{\text{X1}} \approx 0$). Upon lineage commitment in G0 progenitors, cooperative TF regulons are co-activated, causing a significant correlation gain ($\Delta r = r_{\text{G0}} - r_{\text{X1}} > 0$). In empirical single-cell datasets, raw Pearson correlation gains between TF pairs typically span $\Delta r \in [0.10, 0.35]$. Multiplying by $3.0\times$ expands this dynamic range such that a gain of $\Delta r \ge 0.33$ achieves maximum score (1.0), reflecting complete regulatory recruitment.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:654`


### 1.5 Cluster Specificity Formulation ($\text{Inverse Breadth}$)

$$\text{Score}_{\text{specificity}} = \frac{1.0}{n_{\text{clusters supporting}}}$$

- **Biological Rationale**: Genuine neural master regulators are confined to neural and neural-progenitor clusters. Broadly expressed transcriptional machinery, ubiquitous chromatin modifiers, and pleiotropic TFs show multi-cluster expression ($n_{\text{clusters}} \gg 1$). The reciprocal cluster count acts as an explicit penalty against pleiotropic expression.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:421-424`


### 1.6 Perez Lineage TF Superfamily Classification
- **1.0**: Assigned if TF class belongs to Neural Superfamilies (bHLH, Homeobox, POU, C2H2, SOX, FOX, etc.)
- **0.5**: Assigned if TF class is another confirmed structural TF family
- **0.0**: Assigned if unclassified or absent from Perez MOESM5 catalog

- **Biological Rationale**: Perez et al. (*Nat. Commun.* 2025) structurally classified the *S. mediterranea* transcription factor repertoire across differentiation trajectories. Comparative genomics demonstrates that metazoan neural specification is governed by conserved DNA-binding domain families (bHLH, POU-homeodomain, LIM-homeobox, C2H2 zinc fingers). Weighting known neural families at 1.0, general TFs at 0.5, and non-TFs at 0.0 injects established structural biology priors without discarding novel TF families.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:508-550`


### 1.7 Perez Regulatory Influence Score

$$\text{Score}_{\text{perez-influence}} = \text{Influence}_{\text{neuron fate}}$$

- **Biological Rationale**: [Perez et al. (2025)](https://www.nature.com/articles/s41467-025-65712-0#Sec94) computed ANANSE regulatory influence scores across 9 cell fates (MOESM19). The neuron fate influence score represents the normalized rank of each TF's regulatory impact in neural differentiation. A score of 1.0 means the TF has the highest regulatory influence in neuron fate specification.
- **Data Source**: `41467_2025_65712_MOESM19_ESM.xlsx`, sheet `infl_neuron_neoblast_250k` from [Nature Communications](https://www.nature.com/articles/s41467-025-65712-0#Sec94)
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:integrate_perez_influence()`


### 1.8 Evidence Stream Weight Normalization

The 9 evidence streams are weighted as follows:

$$\mathbf{w}_{\text{default}} = [0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]$$

- **Expression** (0.2): Highest weight — direct evidence of transcriptional activity
- **All other 8 streams** (0.1 each): Equal weight — supporting evidence streams

The EvidenceScorer always renormalizes over **present** streams per candidate, so the effective weight depends on which streams have data for each candidate. This ensures that missing evidence does not penalize candidates as *absence of evidence is not evidence of absence*.

### 1.9 Formula Revisions (2026-09-04 mathematical audit)

The following refinements were applied after a full mathematical audit; the
rationales above describe the current (post-audit) behavior:

1. **True log2FC (Fincher/Plass/Cui)** — Scanpy's `rank_genes_groups`
   `logfoldchanges` are a difference of log1p values (pseudo-fold-changes),
   not real log2 fold changes. The pipeline now computes the true log2FC
   from linear-space cluster means (`expm1` of the log1p data, pseudocount
   $10^{-9}$), making the divisor-5.0 cap in §1.1 directly comparable with
   King's mmc7 log2FC.
2. **Best-cluster selection** — the argmax positive fold-change is now taken
   among *significantly* enriched clusters (BH $q \le 0.10$); a gene whose
   largest fold-change cluster is non-significant but which has a
   significant cluster elsewhere is scored on the significant evidence.
3. **King specificity units** — fractional breadth $1-(n-1)/(N-1)$ now uses
   (compartment, subcluster) PAIRS in both numerator and denominator
   (175 distinct subcluster names vs 180 compartment-subcluster pairs; five
   names occur in both G0 and X1).
4. **Perez lineage sentinel** — MOESM5 uses `-` for the ~58k non-TF genes;
   these are now treated as absent (0.0, no atlas membership) instead of
   receiving the 0.5 "other TF class" score.
5. **RNAi/correlation ID matching** — the short-ID parser now extracts the
   numeric gene field from structured IDs (`dd_Smed_v6_11150_0_1 → dd11150`;
   the previous lazy regex returned `dd6`), restoring RNAi-stream matches
   for 18 known neural TFs; correlation Δr uses joint NaN masking so the
6. **King neural gate (2026-09-04 audit)** — the neural-enrichment gate on
   King mmc7 G0 subclusters is `log2FC ≥ 1.5`, matching King 2024's own
   STAR Methods criterion ("at least a 1.5 log2FC enrichment"; every mmc7
   value already passed the authors' upstream p ≤ 0.001 MAST filter). The
   previous stricter gate of 2.0 dropped 7 TFs the paper itself reports as
   neural-enriched (dd17385/sp6-9 1.75, dd19255 1.78, dd36480 1.93,
   dd47123 1.77, dd5882 1.70, dd63520 1.51, dd7442 1.97).
7. **Unused raw atlas arms (2026-09-04 scope decision)** — the following
   downloaded data are deliberately OUT OF SCOPE for the current 5-atlas
   scRNA evidence model and are documented for future extensions:
   - Fincher GSE111764 `BrainClustering` DGE (10,637-cell neuronal
     sub-atlas) and `SexualClustering` DGE — a dedicated Fincher-brain
     evidence stream (neuronal-subtype refinement) is future work.
   - Cui OMIX003867 spatial arm: 12 Visium h5ads (0/6/12 hpa, 1/3/7 dpa,
     OMIX003867-02/-03) plus `adata_Neoblast.h5ad` and `plk1_cut5d.h5ad` —
     spatially-resolved and neoblast-perturbation streams are future work.
    The current pipeline uses 100% of the principal-cells arms: Fincher
    50,562 × 26,561; Plass 37,507 × 28,674; Cui 55,014 × 19,198 (69.4% of
    Cui genes retained — capped by the Rosetta Stone SMED→v6 mapping, not a
    pipeline choice; 7,265 SMED genes have no v6 mapping).
8. **Joint NaN masking (correlation Δr)** — x1/g0 values of a pair row
   can never be re-aligned independently (dropping NaNs per column
   silently re-pairs a different row's values).
9. **HVG forcing** — only the King mmc4 catalog (418 TFs) is forced into
   the highly-variable panel; the 14k master catalog seeds records but
   does not bias clustering.
10. **Perez influence mapping** — restricted to 1:1 reciprocal-best-hit
    v6↔h1SMcG pairs (the collapsed `Similar` column claims 14.4k of 25k v6
    IDs for more than one h1SMcG, making first-wins attribution arbitrary).
11. **Composite score unclipped (2026-09-04 audit)** — composite =
    base + bonuses (max +0.07) with NO 1.0 clip. The former clip
    saturated 6+ genes at exactly 1.0 across methods, destroying ranking
    resolution exactly at the top; deterministic tie-breaks (composite ->
    method base -> integrated -> n_streams -> gene_id) are the sole
    ordering authority.
12. **Consensus null (2026-09-04 audit)** — cross-method consensus uses
    the randomization null p0 = 10/N (N = candidate universe), matching
    overlap_significance.py; the previous binomial p=1/3 null was invalid
    by ~390× and structurally zero-power. Permutation p-values use the
    add-one estimator (b+1)/(n+1); genes whose scores are fully explained
    by the label-independent King table are flagged
    `untestable_by_permutation` (their cluster-label permutation p is
    1.0 by construction, and the honest statement is that their rank
    rests on external tables, not cluster-specific expression).

---


## 2. Statistical Testing & False Discovery Rate Control

### 2.1 Benjamini-Hochberg Multiple Testing Correction ($\text{FDR } q \le 0.10$)

- **Statistical Rationale**: In differential expression analysis across single-cell clusters (Wilcoxon rank-sum test), thousands of simultaneous gene-cluster hypotheses are evaluated. We apply the Benjamini-Hochberg step-up procedure to control the False Discovery Rate (FDR). An exploratory threshold of $q \le 0.10$ is chosen to avoid type II errors (false negatives) during initial candidate harvesting, while downstream multi-stream integration eliminates spurious single-test hits.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:361`

---

## 3. Bayesian Uncertainty Quantification (Dirichlet Sampling)

### 3.1 Centered Dirichlet Concentration Parameter ($k = 40.0$)

$$\mathbf{w}^{(m)} \sim \text{Dirichlet}(k \cdot \mathbf{w}_{\text{default}}), \quad m = 1, \dots, 1000$$

- **Mathematical Rationale**: Fixed-weight scoring models assume exact certainty in parameter weights. The centered Dirichlet model formalizes weight uncertainty by drawing 1,000 weight vectors centered on $\mathbf{w}_{\text{default}} = [0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]$. Setting the concentration parameter to $k = 40.0$ corresponds to 40 pseudo-observations of evidence reliability, yielding a 95% credible interval of approximately $\pm 0.10$ around each baseline weight.
- **Location**: `projects/NeuralTF/scripts/dirichlet_centered.py` (evaluated across all candidates in `rank.csv`)

### 3.2 Uniform Dirichlet Prior ($\alpha_i = 1.0, \, \forall i$)

$$\mathbf{w}^{(m)} \sim \text{Dirichlet}(\mathbf{1}_9)$$

- **Mathematical Rationale**: To demonstrate that candidate rankings are driven by intrinsic biological signal rather than investigator weight choices, the uniform Dirichlet samples uniformly across the 9-dimensional probability simplex ($\alpha_i = 1$). Concordance between uniform Dirichlet medians and fixed-weight rankings confirms high stability across all candidates without prior weighting assumptions.
- **Location**: `projects/NeuralTF/scripts/dirichlet_uniform.py` (evaluated across all candidates in `rank.csv`)

---

## 4. Candidate Selection Funnel

```
┌─────────────────────────────────────────────────────────────────┐
│ King 2024 Master TF Catalog (mmc4.xlsx)                         │
│ 418 Candidate TFs seeded (TF? != NA)                            │
└───────────────────────────────┬─────────────────────────────────┘
                                │ Wilcoxon DE across Leiden clusters
                                │ in Fincher, Plass, and Cui atlases
                                ▼ (FDR q ≤ 0.10)
┌─────────────────────────────────────────────────────────────────┐
│ TFs with significant single-cell cluster DE                     │
└───────────────────────────────┬─────────────────────────────────┘
                                 │ Integration of King 2024 mmc7
                                 │ G0 progenitor neural subclusters (log2FC ≥ 1.5)
                                 │ + King mmc5 RNAi screen targets
                                │ + King mmc6 TF-pair correlations
                                ▼
┌─────────────────────────────────────────────────────────────────┐
│ TFs scored across King evidence streams                         │
└───────────────────────────────┬─────────────────────────────────┘
                                │ Integration of Perez 2025:
                                │ - MOESM5 TF lineage classification
                                │ - MOESM19 ANANSE regulatory influence
                                ▼
┌─────────────────────────────────────────────────────────────────┐
│ 278 Total Candidate TFs (rank.csv)                              │
│ Scored across all 9 evidence streams                            │
└───────────────────────────────┬─────────────────────────────────┘
                                │ Neural Gate:
                                │ (neural_enriched > 0) | (rnai > 0)
                                ▼
┌─────────────────────────────────────────────────────────────────┐
│ 101 Neural-Enriched Candidate TFs (rank_neural.csv)             │
│ (Dual-track shortlist: Track A + Track B)                       │
└─────────────────────────────────────────────────────────────────┘
```
