# NeuralTF Pipeline — Parameter Scientific Rationale & Sensitivity Documentation

This document records the biological, statistical, and mathematical justifications for all parameters, thresholds, and scoring formulations used in the `NeuralTF` pipeline, formatted according to *Nature Communications* reporting standards.

---

## 1. Multi-Atlas Evidence Stream Scoring Formulations

### 1.1 Expression Score Cap ($\text{Divisor} = 5.0$)

$$\text{Score}_{\text{expression}} = \min\left(1.0, \frac{\max(\text{log}_2\text{FC})}{5.0}\right)$$

- **Biological Rationale**: In planarian single-cell RNA sequencing (*Schmidtea mediterranea*), a $\text{log}_2\text{FC} = 5.0$ corresponds to a $2^5 = 32$-fold linear upregulation relative to background non-target cells. Transcription factor (TF) steady-state transcript abundance in whole-animal or tissue clusters reaches saturation in target promoter occupancy well before 32-fold upregulation. Setting the divisor at 5.0 maps the physiological range $[0, 5.0]$ linearly to $[0, 1.0]$. This avoids outlier compression where a single extreme fold-change (e.g., $\text{log}_2\text{FC} = 12$) would compress biologically meaningful, moderately expressed lineage master regulators ($\text{log}_2\text{FC} \in [2.0, 4.0]$) into near-zero scores.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:413-416`
- **Sensitivity Analysis**:
  | Divisor | Total Candidates | Top-10 Jaccard vs 5.0 | Biological Behavior |
  |:---:|:---:|:---:|:---|
  | 3.0 | 289 | 0.88 | Premature saturation at 8-fold FC; compresses high-expression dynamic range |
  | **5.0** | **289** | **1.00** | **Optimal linear dynamic range (Default)** |
  | 8.0 | 289 | 0.95 | Down-weights moderately expressed but functional TFs |
  | 10.0 | 289 | 0.90 | Biases towards high-copy transcripts at the expense of regulatory factors |

---

### 1.2 Multi-Atlas Reproducibility Denominator ($N_{\text{atlases}} = 5$)

$$\text{Score}_{\text{reproducibility}} = \frac{\min(n_{\text{supporting atlases}}, 5)}{5.0}$$

- **Biological Rationale**: True lineage specification drivers exhibit consistent transcriptional induction across independent experimental protocols, cell dissociation methods, and sequencing technologies. The pipeline interrogates all 5 single-cell/regulatory atlases:
  1. **Fincher 2018** (Drop-seq, whole-animal, dd_Smed_v4 transcript model)
  2. **Plass 2018** (Drop-seq, whole-animal, dd_Smed_v6 transcript model)
  3. **Cui 2023** (10x Genomics Chromium, regeneration single-cell atlas, 55,014 cells)
  4. **King 2024** (FACS-purified G0/X1 progenitor single-cell atlas & TF catalog)
  5. **Perez 2025** (Lineage Single-Cell Differentiation Atlas & TF Domain Family Classification)
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:663-672`
- **Scientific Impact**: Candidate TFs supported across all 5 atlases achieve $s_{\text{repro}} = 1.0$, while single-atlas detections receive $s_{\text{repro}} = 0.20$.

---

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

---

### 1.4 Co-Expression Correlation Gain Multiplier ($\text{Multiplier} = 3.0$)

$$\text{Score}_{\text{correlation}} = \min\left(1.0, \max(0.0, r_{\text{G0}} - r_{\text{X1}}) \times 3.0\right)$$

- **Biological Rationale**: In pluripotent X1 neoblasts, neural TF pairs display negligible co-expression correlation ($r_{\text{X1}} \approx 0$). Upon lineage commitment in G0 progenitors, cooperative TF regulons are co-activated, causing a significant correlation gain ($\Delta r = r_{\text{G0}} - r_{\text{X1}} > 0$). In empirical single-cell datasets, raw Pearson correlation gains between TF pairs typically span $\Delta r \in [0.10, 0.35]$. Multiplying by $3.0\times$ expands this dynamic range such that a gain of $\Delta r \ge 0.33$ achieves maximum score (1.0), reflecting complete regulatory recruitment.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:654`

---

### 1.5 Cluster Specificity Formulation ($\text{Inverse Breadth}$)

$$\text{Score}_{\text{specificity}} = \frac{1.0}{n_{\text{clusters supporting}}}$$

- **Biological Rationale**: Genuine neural master regulators are confined to neural and neural-progenitor clusters. Broadly expressed transcriptional machinery, ubiquitous chromatin modifiers, and pleiotropic TFs show multi-cluster expression ($n_{\text{clusters}} \gg 1$). The reciprocal cluster count acts as an explicit penalty against pleiotropic expression.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:421-424`

---

### 1.6 Perez Lineage TF Superfamily Classification

$$\text{Score}_{\text{perez\_lineage}} = \begin{cases} 1.0 & \text{if TF class } \in \text{Neural Superfamilies (bHLH, Homeobox, POU, C2H2, SOX, FOX, etc.)} \\ 0.5 & \text{if TF class is another confirmed structural TF family} \\ 0.0 & \text{if unclassified or absent from Perez MOESM5 catalog} \end{cases}$$

- **Biological Rationale**: Perez et al. (*Nat. Commun.* 2025) structurally classified the *S. mediterranea* transcription factor repertoire across differentiation trajectories. Comparative genomics demonstrates that metazoan neural specification is governed by conserved DNA-binding domain families (bHLH, POU-homeodomain, LIM-homeobox, C2H2 zinc fingers). Weighting known neural families at 1.0, general TFs at 0.5, and non-TFs at 0.0 injects established structural biology priors without discarding novel TF families.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:508-550`

---

## 2. Statistical Testing & False Discovery Rate Control

### 2.1 Benjamini-Hochberg Multiple Testing Correction ($\text{FDR } q \le 0.10$)

- **Statistical Rationale**: In differential expression analysis across single-cell clusters (Wilcoxon rank-sum test), thousands of simultaneous gene-cluster hypotheses are evaluated. We apply the Benjamini-Hochberg step-up procedure to control the False Discovery Rate (FDR). An exploratory threshold of $q \le 0.10$ is chosen to avoid type II errors (false negatives) during initial candidate harvesting, while downstream multi-stream integration eliminates spurious single-test hits.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:361`

---

## 3. Bayesian Uncertainty Quantification (Dirichlet Sampling)

### 3.1 Centered Dirichlet Concentration Parameter ($k = 40.0$)

$$\mathbf{w}^{(m)} \sim \text{Dirichlet}(k \cdot \mathbf{w}_{\text{default}}), \quad m = 1, \dots, 1000$$

- **Mathematical Rationale**: Fixed-weight scoring models assume exact certainty in parameter weights. The centered Dirichlet model formalizes weight uncertainty by drawing 1,000 weight vectors centered on $\mathbf{w}_{\text{default}} = [0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]$. Setting the concentration parameter to $k = 40.0$ corresponds to 40 pseudo-observations of evidence reliability, yielding a 95% credible interval of approximately $\pm 0.10$ around each baseline weight.
- **Location**: `projects/NeuralTF/scripts/dirichlet_prioritize.py:71`, `dirichlet_centered_all249.py`

### 3.2 Uniform Dirichlet Prior ($\alpha_i = 1.0, \, \forall i$)

$$\mathbf{w}^{(m)} \sim \text{Dirichlet}(\mathbf{1}_8)$$

- **Mathematical Rationale**: To prove that candidate rankings are driven by intrinsic biological signal rather than investigator weight choices, the uniform Dirichlet samples uniformly across the 8-dimensional probability simplex ($\alpha_i = 1$). Concordance between uniform Dirichlet medians and fixed-weight rankings (**10/10 top-10 overlap**) confirms extreme robustness.
- **Location**: `projects/NeuralTF/scripts/dirichlet_uniform.py`, `dirichlet_uniform_all249.py`

---

## 4. Candidate Selection Funnel Audit

```
┌─────────────────────────────────────────────────────────────────┐
│ King 2024 Master TF Catalog (mmc4.xlsx)                         │
│ 418 Candidate TFs seeded (TF? != NA)                            │
└───────────────────────────────┬─────────────────────────────────┘
                                │ Wilcoxon DE across Leiden clusters
                                │ in Fincher, Plass, and Cui atlases
                                ▼ (FDR q ≤ 0.10)
┌─────────────────────────────────────────────────────────────────┐
│ 224 TFs with significant single-cell cluster DE                 │
└───────────────────────────────┬─────────────────────────────────┘
                                │ Integration of King 2024 mmc7
                                │ G0 progenitor neural subclusters (log2FC ≥ 2)
                                ▼ + King mmc5 RNAi screen targets
┌─────────────────────────────────────────────────────────────────┐
│ 289 Total Candidate TFs (rank.csv)                              │
│ Scored across all 8 evidence streams                            │
└───────────────────────────────┬─────────────────────────────────┘
                                │ Neural Gate:
                                │ (neural_enriched > 0) | (rnai > 0)
                                ▼
┌─────────────────────────────────────────────────────────────────┐
│ 102 Neural-Enriched Candidate TFs (rank_neural.csv)             │
│ (96 G0 neural subcluster hits ∪ 6 RNAi-validated TFs)           │
└─────────────────────────────────────────────────────────────────┘
```