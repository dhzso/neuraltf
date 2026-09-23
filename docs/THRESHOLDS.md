# NeuralTF Pipeline — Parameter Scientific Rationale & Sensitivity Documentation

This document records the biological, statistical, and mathematical justifications for all parameters, thresholds, and scoring formulations used in the `NeuralTF` pipeline.

> **Alignment note (2026-09-23)**: every formula and code reference below was re-verified against the checked-in implementation (`src/bioforge/projects/neuraltf/pipeline.py`, `src/bioforge/evidence/scoring.py`) and artifacts (`projects/NeuralTF/runs/pipeline_run/rank.csv`, `projects/NeuralTF/data/king_atlas.tsv`); the corrections made in this pass are recorded in §1.9 item 13.

---

## 1. Multi-Atlas Evidence Stream Scoring Formulations

### 1.0 Symbol Reference (how to read the formulas)

| Symbol | Meaning | Section |
|:---|:---|:---|
| $n_{\text{sig clusters}}$ | number of Leiden clusters of one atlas in which the gene is *significantly* upregulated (BH $q \le 0.10$ and positive fold-change) | §1.5 (atlas component) |
| $n_{\text{King pairs}}$ | number of distinct (compartment, subcluster) pairs with a King-atlas hit for the gene | §1.5 (King component) |
| $N_{\text{King}} = 181$ | all King (compartment, subcluster) pairs (176 distinct subcluster names) | §1.5, §1.9 item 3 |
| $n_{\text{neural pairs}}$, $N_{\text{neural}} = 80$ | the same pair count restricted to the neural mask | §1.5a |
| $k = 40$ | Dirichlet concentration (pseudo-observation strength) — **unrelated to clusters** | §3.1 |
| $k_{\text{clusters}}$ | number of clusters a gene was *tested* in; used only as the Sidak exponent $1-(1-p)^{k}$ when un-conditioning stored p-values | `scripts/stats/meta_analytic_pvalue.py` |

**"Breadth"** means how many clusters/pairs the TF is (significantly) expressed in. Both specificity
formulations invert it, but with different shapes: the principal-atlas component uses the reciprocal
$1/n$ (a hard cliff: $n=1 \to 1.0$, $n=2 \to 0.5$), while the King/neural components use the linear
rescaling $1-(n-1)/(N-1)$ ($n=1 \to 1.0$, $n=181 \to 0.0$; smooth, no cliff). The recorded stream score
is the **maximum** over the components (§1.5), and evidence-card notes record only the **last** component
written — so a card note can disagree with the value (e.g. `specificity: 1.00 — atlas=cui,n_sig_clusters=9`
means another atlas already scored 1.0 and the note was then overwritten, not that cui contributed 1.0).

### 1.1 Expression Score Cap ($\text{Divisor} = 5.0$)

$$\text{Score}_{\text{expression}} = \min\left(1.0, \frac{\max(\text{log}_2\text{FC})}{5.0}\right)$$

- **Biological Rationale**: In planarian (*Schmidtea mediterranea*) single-cell RNA (scRNA) sequencing, a $\text{log}_2\text{FC} = 5.0$ corresponds to a $2^5 = 32$-fold of linear upregulation which is relative to background non-target cells. Sometimes, transcription factor (TF) steady-state transcript abundance in whole-animal or tissue clusters reaches saturation in target promoter occupancy well before 32-fold upregulation. Setting the divisor at 5.0 maps the physiological range $[0, 5.0]$ linearly to $[0, 1.0]$. This avoids outlier compression where a single extreme fold-change (e.g., $\text{log}_2\text{FC} = 12$) would compress biologically meaningful, moderately expressed lineage master regulators ($\text{log}_2\text{FC} \in [2.0, 4.0]$) into near-zero scores.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:793-802` (Fincher/Plass/Cui), `pipeline.py:851-854` + `pipeline.py:900` (King)
- **Sensitivity Analysis**: the divisor rescales only the expression stream, so it cannot change candidate *eligibility* — every record is written to `rank.csv` (11,696 candidates) and the neural gate (`neural_enriched > 0 | rnai > 0`) never reads the expression score. The divisor controls only where the expression stream saturates:
  | Divisor | Saturation Point | Effect on the Expression Stream |
  |:---:|:---:|:---|
  | 3.0 | $2^3 = 8$-fold FC | Saturation arrives early; compresses the upper dynamic range |
  | **5.0** | **$2^5 = 32$-fold FC (Default)** | **Maps the physiological range $[0, 5]$ linearly to $[0, 1]$** |
  | 8.0 | $2^8 = 256$-fold FC | Down-weights moderately expressed but functional TFs |
  | 10.0 | $2^{10} = 1024$-fold FC | Biases towards high-copy transcripts at the expense of regulatory factors |

  (An earlier revision of this table reported 289 total candidates and top-10 Jaccard values; neither is reproducible from the current 5-atlas / 11-stream configuration — the candidate universe does not depend on the divisor.)


### 1.2 Multi-Atlas Reproducibility Denominator ($N_{\text{atlases}} = 5$)

$$\text{Score}_{\text{reproducibility}} = \frac{\min(n_{\text{supporting atlases}}, 5)}{5.0}$$

- **Biological Rationale**: True lineage specification drivers exhibit consistent transcriptional induction across independent experimental protocols, cell dissociation methods, and sequencing technologies. The pipeline interrogates all 5 single-cell/regulatory atlases:
  1. **Fincher 2018** (Drop-seq, whole-animal, `dd_Smed_v4` transcript model)
  2. **Plass 2018** (Drop-seq, whole-animal, `dd_Smed_v6` transcript model)
  3. **Cui 2023** (10x Genomics Chromium, regeneration single-cell atlas, 55,014 cells)
  4. **King 2024** (FACS-purified G0/X1 progenitor single-cell atlas & TF catalog)
  5. **Perez 2025** (Lineage Single-Cell Differentiation Atlas & TF Domain Family Classification)
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:1276-1286` (`assign_reproducibility`)
- **Scientific Impact**: Candidate TFs supported across all 5 atlases achieve $s_{\text{repro}} = 1.0$, while single-atlas detections receive $s_{\text{repro}} = 0.20$.


### 1.3 King Atlas Neural Subcluster Enrichment Threshold ($\text{log}_2\text{FC} \ge 1.5$)

$$\text{Neural Gate} = \mathbb{I}\left(\text{subcluster} \in \text{Neural} \land \text{log}_2\text{FC} \ge 1.5\right)$$

- **Biological Rationale**: In planarian adult stem cell differentiation (King et al., *Cell Reports* 2024), pluripotent neoblasts (X1) exit the cell cycle into post-mitotic committed progenitors (G0). Fate-specifying transcription factors (FSTFs) undergo sharp, switch-like transcriptional activation during this transition. A $\text{log}_2\text{FC} \ge 1.5$ represents a $\ge 2.83$-fold enrichment in neural G0 progenitors compared to general G0 populations, matching the King 2024 STAR Methods criterion ("at least a 1.5 log2FC enrichment"); every mmc7 value already passed the authors' upstream $p \le 0.001$ MAST filter. The previous stricter gate of 2.0 dropped 8 TFs the paper itself reports as neural-enriched (see §1.9 item 6).
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:47` (`_NEURAL_FC_THRESHOLD = 1.5`)
- **Sensitivity Analysis**:
  | $\text{log}_2\text{FC}$ Cutoff | King Neural Genes | `rank_neural` Gate (∪ RNAi) | Neural-Set Jaccard vs 1.5 | Statistical Assessment |
  |:---:|:---:|:---:|:---:|:---|
  | 1.0 | 116 | 145 | 0.98 | Permissive; admits sub-threshold progenitors the paper does not call neural-enriched |
  | **1.5** | **114** | **143** | **1.00** | **Matches King 2024 STAR Methods criterion (Default)** |
  | 2.0 | 106 | 135 | 0.93 | Stricter than the paper; drops the 8 TFs in the $[1.5, 2.0)$ band |
  | 2.5 | 98 | 130 | 0.86 | Overly stringent; excludes validated neural TFs with modest basal expression |
  | 3.0 | 88 | 122 | 0.77 | Severe loss of known RNAi-validated regulators |

  Counts recomputed (2026-09-23) from `data/king_atlas.tsv` — "King Neural Genes" = unique genes whose maximum log2FC within the neural mask (G0 `neural*`, X1 `cell_type == Neural`, X1 Major Tissue `Neural`) clears the cutoff — and `runs/pipeline_run/rank.csv` (the 68 RNAi-screened genes). "Neural-Set Jaccard" compares the unique King neural gene set at each cutoff against the 1.5 default.


### 1.4 Co-Expression Correlation Gain Multiplier ($\text{Multiplier} = 3.0$)

$$\text{Score}_{\text{correlation}} = \min\left(1.0, \max(0.0, r_{\text{G0}} - r_{\text{X1}}) \times 3.0\right)$$

- **Biological Rationale**: In pluripotent X1 neoblasts, neural TF pairs display negligible co-expression correlation ($r_{\text{X1}} \approx 0$). Upon lineage commitment in G0 progenitors, cooperative TF regulons are co-activated, causing a significant correlation gain ($\Delta r = r_{\text{G0}} - r_{\text{X1}} > 0$). In empirical single-cell datasets, raw Pearson correlation gains between TF pairs typically span $\Delta r \in [0.10, 0.35]$. Multiplying by $3.0\times$ expands this dynamic range such that a gain of $\Delta r \ge 0.33$ achieves maximum score (1.0), reflecting complete regulatory recruitment.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:1253-1268` (`integrate_correlations`; the ×3.0 gain multiplier is applied at line 1264)


### 1.5 Cluster Specificity Formulations ($\text{Inverse Breadth}$)

$$\text{Score}_{\text{specificity}} = \max\left(\underbrace{\frac{1.0}{n_{\text{sig clusters}}}}_{\text{principal atlases}},\; \underbrace{1 - \frac{n_{\text{King pairs}} - 1}{181 - 1}}_{\text{King atlas}}\right)$$

- **Principal atlases (Fincher/Plass/Cui)** — computed per atlas: inverse of the number of *significantly* upregulated Leiden clusters (BH $q \le 0.10$ and positive fold-change) in that atlas, so the atlas with the narrowest significant footprint sets this component (a TF significant in 1 cluster scores 1.0; one significant in every cluster scores $1/n_{\text{clusters}}$) (`pipeline.py:804-813`).
- **King atlas** — fractional breadth $1-(n-1)/(N-1)$ over (compartment, subcluster) pairs with $N = 181$ pairs (same units as §1.9 item 3); a smooth, graded penalty that avoids the steep $1 \to 2$-cluster cliff of $1/n$ (`pipeline.py:901-911`).
- The recorded stream score is the **maximum** of the atlas and King components (each `add_score` call passes `max(existing, new)`).
- **Biological Rationale**: Genuine neural master regulators are confined to neural and neural-progenitor clusters. Broadly expressed transcriptional machinery, ubiquitous chromatin modifiers, and pleiotropic TFs show multi-cluster expression ($n_{\text{clusters}} \gg 1$). The reciprocal cluster count (and the King fractional breadth) act as an explicit penalty against pleiotropic expression.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:804-813` (atlas $1/n$), `pipeline.py:901-911` (King fractional breadth)

#### 1.5a Neural Specificity (King neural subclusters)

$$\text{Score}_{\text{neural-specificity}} = 1 - \frac{n_{\text{neural pairs}} - 1}{80 - 1}$$

- Fractional breadth over the (compartment, subcluster) pairs that belong to the neural mask; $N_{\text{neural}} = 80$ pairs. Scored only for genes with a neural-enriched hit at the $\text{log}_2\text{FC} \ge 1.5$ gate (§1.3): a gene hit in a single neural pair scores 1.0, one spread across all 80 scores 0.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:919-926`


### 1.6 Perez Lineage TF Superfamily Classification
- **1.0**: Assigned if the Perez `tf_class` matches (case-insensitive substring) a neural-relevant family — the exact `_PEREZ_NEURAL_CLASSES` set: bHLH, Homeodomain, LHX, POU, Forkhead, PAX, NKX, SIX, EGR, GLI, IRX, COE, HMG (SOX-class), HOX, ISL.
- **0.5**: Assigned if TF class is another confirmed structural TF family (e.g. generic `zf-C2H2`, Tbox, SMAD, GCM, PBX — **not** neural-weighted)
- **0.0**: Assigned if unclassified or absent from Perez MOESM5 catalog

- **Biological Rationale**: Perez et al. (*Nat. Commun.* 2025) structurally classified the *S. mediterranea* transcription factor repertoire across differentiation trajectories. Comparative genomics demonstrates that metazoan neural specification is governed by conserved DNA-binding domain families (bHLH, POU-homeodomain, LIM-homeobox (LHX/ISL), HMG/SOX, Forkhead, PAX, NKX, SIX, EGR, GLI, IRX, COE, HOX). Weighting known neural families at 1.0, general TFs at 0.5, and non-TFs at 0.0 injects established structural biology priors without discarding novel TF families; a confirmed class (score ≥ 0.5) — and only a confirmed class — also counts as Perez atlas support for the reproducibility denominator.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:962-978` (class set), `pipeline.py:980-1024` (`integrate_perez`)


### 1.7 Perez Regulatory Influence Score

$$\text{Score}_{\text{perez-influence}} = \text{Influence}_{\text{neuron fate}}$$

- **Biological Rationale**: [Perez et al. (2025)](https://www.nature.com/articles/s41467-025-65712-0#Sec94) computed ANANSE regulatory influence scores across 9 cell fates (MOESM19). The neuron fate influence score represents the normalized rank of each TF's regulatory impact in neural differentiation. A score of 1.0 means the TF has the highest regulatory influence in neuron fate specification.
- **Data Source**: `41467_2025_65712_MOESM19_ESM.xlsx`, sheet `infl_neuron_neoblast_250k` from [Nature Communications](https://www.nature.com/articles/s41467-025-65712-0#Sec94)
- **What is actually scored**: the ANANSE *influence score* of the TF on the **neuron fate** — how strongly the TF's inferred enhancer-mediated regulatory network explains the neuron (vs neoblast) expression program. In that sheet the `influence_score` column is a **rank-scaled** quantity: across its 92 factors it equals exactly $1 - i/(n-1)$ (0-based rank $i$, $n = 92$, step $1/91 \approx 0.011$), so 1.0 = top-ranked neuron-fate regulator and 0.0 = lowest. The pipeline uses `influence_score` only; the raw network quantity (`influence_score_raw`) and its components (`target_score`, `G_score`, `direct_targets`, `factor_fc`) are not used.
- **Coverage**: sheet factors are h1SMcG (h1 genome) IDs and are mapped to v6 by 1:1 reciprocal-best-hit pairs only (§1.9 item 10), so the stream is present on just **71 of 11,696** candidates in the current run (70 with a score > 0). A candidate without the stream is **not** penalised — the scorer renormalizes over present streams (§1.8).
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:1030-1111` (`integrate_perez_influence`; the sheet is auto-discovered as the first MOESM19 sheet containing `neuron`)


### 1.8 Evidence Stream Weight Normalization

The 11 evidence streams are weighted as follows (single source of truth:
`bioforge.evidence.scoring.DEFAULT_WEIGHTS` / `STREAM_ORDER`, declared in
`src/bioforge/evidence/scoring.py:26-60`):

$$\mathbf{w}_{\text{default}} = [0.10, 0.10, 0.10, 0.05, 0.05, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10]$$

- **Expression, Specificity, Reproducibility, Neural Enriched, Neural Specificity,
  Perez Lineage, Perez Influence, Fincher Brain, Cui Temporal** (0.10 each) — the
  independent or widely-supported evidence streams.
- **RNAi, Correlation** (0.05 each) — the two streams that encode the ground-truth
  label (King mmc5/mmc6); deliberately de-emphasized so no single label stream
  saturates the integrated score.

The EvidenceScorer always renormalizes over **present** streams per candidate, so the effective weight depends on which streams have data for each candidate. This ensures that missing evidence does not penalize candidates as *absence of evidence is not evidence of absence*.

#### 1.8a Fincher Brain stream (`fincher_brain`, w=0.10)

$$\text{Score}_{\text{fincher-brain}} = \min\left(1.0, \frac{\max(\text{log}_2\text{FC}_{\text{brain clusters}})}{5.0}\right)$$

Computed from the **independent** Fincher 2018 `BrainClustering` DGE (7,766 head cells,
26,565 v4 genes; `scripts/convert_fincher_brain.py`), de-novo Leiden clustered and
scored exactly like the principal expression stream (Wilcoxon DE, global BH-FDR
q ≤ 0.10, true log2FC via `_cluster_log2fc`). It does **not** count toward the
reproducibility denominator (which stays at the 5 original atlases).

#### 1.8b Cui Temporal stream (`cui_temporal`, w=0.10)

$$\text{Score}_{\text{cui-temporal}} = \min\left(1.0, \frac{\log_2\!\left(\frac{\text{peak} + 0.1}{\text{baseline} + 0.1}\right)}{2.0}\right)$$

Computed in `BigCellType == "Neuronal"` cells across the 8 Cui 2023 regeneration
timepoints (cut0d ... cut7d), as the log2 fold-change of the peak post-amputation
mean over the cut0d baseline. Genes below a 0.1 (linear counts-per-10k) neuronal
expression floor receive score 0 (no temporal claim from near-zero expression).

### 1.9 Formula Revisions

The following refinements define the current behavior:

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
3. **King specificity units** — fractional breadth $1-(n-1)/(N-1)$ uses
   (compartment, subcluster) PAIRS in both numerator and denominator:
   181 compartment-subcluster pairs (127 G0 Progenitor + 46 X1 + 8 X1
   Major Tissue) collapse to 176 distinct subcluster names, and five names
   occur in more than one compartment. The neural-specificity denominator
   is the 80 neural pairs (§1.5a).
4. **Perez lineage sentinel** — MOESM5 uses `-` for the ~58k non-TF genes;
   these are now treated as absent (0.0, no atlas membership) instead of
   receiving the 0.5 "other TF class" score.
5. **RNAi/correlation ID matching** — the short-ID parser now extracts the
   numeric gene field from structured IDs (`dd_Smed_v6_11150_0_1 → dd11150`;
   the previous lazy regex returned `dd6`), restoring RNAi-stream matches
   for 18 known neural TFs.
6. **King neural gate** — the neural-enrichment gate on the King mmc7
   neural mask (§1.3: G0 `neural*`, X1 `cell_type == Neural`, X1 Major
   Tissue `Neural`) is `log2FC ≥ 1.5`, matching King 2024's own STAR
   Methods criterion ("at least a 1.5 log2FC enrichment"; every mmc7
   value already passed the authors' upstream p ≤ 0.001 MAST filter).
   In the current checked-in atlas, 114 genes clear 1.5; the previous
   stricter gate of 2.0 would drop the 8 TFs in the [1.5, 2.0) band
   (dd7033 1.93, dd14937 1.87, dd9287 1.85, dd11930 1.76, dd5828 1.70,
   dd7442 1.68, dd7227 1.67, dd15473 1.61).
7. **Raw atlas arms (scope decision)** — the following downloaded data are
   now IN SCOPE as dedicated streams (§1.8a/1.8b): Fincher `BrainClustering`
   DGE (`fincher_brain`), and Cui regeneration time-course (`cui_temporal`).
   Remaining FUTURE work (documented, not integrated):
   - Fincher `SexualClustering` DGE.
   - Cui OMIX003867 spatial arm: 12 Visium spatial objects (6 Scanpy
     `.h5ad` in OMIX003867-03: 0/6/12 hpa, 1/3/7 dpa; 6 Seurat `.Robj` in
     OMIX003867-02: cut0h–cut7d) plus `adata_Neoblast.h5ad` and
     `plk1_cut5d.h5ad` — spatially-resolved and neoblast-perturbation
     streams are future work.
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
11. **Composite score unclipped** — composite =
    base + bonuses (max +0.07) with NO 1.0 clip. The former clip
    saturated 6+ genes at exactly 1.0 across methods, destroying ranking
    resolution exactly at the top; deterministic tie-breaks (composite ->
    method base -> integrated -> n_streams -> gene_id) are the sole
    ordering authority.
12. **Consensus null (stratified)** — cross-method consensus tests
    (`overlap_significance.py`, `cross_method_correction.py`) use the
    dual-track stratified null: each method draws 5 genes within the
    tested stratum and 5 within the not-tested stratum, so the per-gene
    null is $5/N_{\text{stratum}}$ (currently 5/68 = 0.0735 tested,
    5/11,628 = 4.3e-4 not-tested; expected pairwise overlap ≈ 0.370)
    instead of the uniform 10/N = 8.5e-4. The legacy uniform-subset
    hypergeometric p (10/N) is retained only as an upper bound — it
    overstates significance by ~25 orders of magnitude for these
    stratified shortlists, and the older binomial p = 1/3 null was
    invalid by ~390× by comparison. Permutation p-values use the
    add-one estimator (b+1)/(n+1); genes whose scores are fully explained
    by the label-independent King table are flagged
    `untestable_by_permutation` (their cluster-label permutation p is
    1.0 by construction, and the honest statement is that their rank
    rests on external tables, not cluster-specific expression).
13. **Documentation alignment (2026-09-23)** — this file was re-synced
    against the checked-in code and artifacts: §1.5 now states the
    actual specificity rule (maximum of the atlas $1/n_{\text{sig}}$ and
    the King fractional breadth) and §1.5a documents the
    neural-specificity formula; the King denominators are 181 (all) and
    80 (neural) pairs over 176 distinct subcluster names; the §1.1 and
    §1.3 sensitivity tables were replaced with values reproducible from
    `rank.csv` / `king_atlas.tsv` (the previous 289-candidate table
    predated the 5-atlas / 11-stream configuration); the Perez
    neural-family list in §1.6 is the exact `_PEREZ_NEURAL_CLASSES` set;
    §1.9 item 12 describes the stratified consensus null; §1.0 adds a
    symbol reference and §1.7 states the rank-scaling and coverage of the
    Perez influence stream; and all code line references were refreshed.

---


## 2. Statistical Testing & False Discovery Rate Control

### 2.1 Benjamini-Hochberg Multiple Testing Correction ($\text{FDR } q \le 0.10$)

- **Statistical Rationale**: In differential expression analysis across single-cell clusters (Wilcoxon rank-sum test), thousands of simultaneous gene-cluster hypotheses are evaluated. We apply the Benjamini-Hochberg step-up procedure to control the False Discovery Rate (FDR). An exploratory threshold of $q \le 0.10$ is chosen to avoid type II errors (false negatives) during initial candidate harvesting, while downstream multi-stream integration eliminates spurious single-test hits.
- **Location**: `src/bioforge/projects/neuraltf/pipeline.py:51` (`_FDR_THRESHOLD = 0.1`), applied at `pipeline.py:523` (Fincher brain), `pipeline.py:708` (principal atlases), and `pipeline.py:745` (King)

---

## 3. Bayesian Uncertainty Quantification (Dirichlet Sampling)

### 3.1 Centered Dirichlet Concentration Parameter ($k = 40.0$)

$$\mathbf{w}^{(m)} \sim \text{Dirichlet}(k \cdot \mathbf{w}_{\text{default}}), \quad m = 1, \dots, 1000$$

- **Mathematical Rationale**: Fixed-weight scoring models assume exact certainty in parameter weights. The centered Dirichlet model formalizes weight uncertainty by drawing 1,000 weight vectors centered on $\mathbf{w}_{\text{default}} = [0.100, 0.100, 0.100, 0.050, 0.050, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100]$. Setting the concentration parameter to $k = 40.0$ corresponds to 40 pseudo-observations of evidence reliability, yielding 95% marginal credible intervals of $[0.029, 0.209]$ for the nine 0.100-weight streams and $[0.006, 0.135]$ for the two 0.050-weight streams (marginal Beta quantiles).
- **Location**: `projects/NeuralTF/scripts/dirichlet_centered.py` (evaluated across all candidates in `rank.csv`; constants `K_DIR = 40.0`, `N_DRAWS = 1000`, `SEED = 2024`, weights imported from `DEFAULT_WEIGHTS`)

### 3.2 Uniform Dirichlet Prior ($\alpha_i = 1.0, \, \forall i$)

$$\mathbf{w}^{(m)} \sim \text{Dirichlet}(\mathbf{1}_{11})$$

- **Mathematical Rationale**: To demonstrate that candidate rankings are driven by intrinsic biological signal rather than investigator weight choices, the uniform Dirichlet samples uniformly across the 11-dimensional probability simplex ($\alpha_i = 1$). Concordance between uniform Dirichlet medians and fixed-weight rankings confirms high stability across all candidates without prior weighting assumptions.
- **Location**: `projects/NeuralTF/scripts/dirichlet_uniform.py` (evaluated across all candidates in `rank.csv`; `alpha = np.ones(n_streams)`, 1,000 draws, seed 2024)

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
│ 11,696 Total Candidate TFs (rank.csv)                              │
│ Scored across all 11 evidence streams                            │
└───────────────────────────────┬─────────────────────────────────┘
                                │ Neural Gate:
                                │ (neural_enriched > 0) | (rnai > 0)
                                ▼
┌─────────────────────────────────────────────────────────────────┐
│ 143 Neural-Enriched Candidate TFs (rank_neural.csv)             │
│ (Dual-track shortlist: tested + not-tested)                       │
└─────────────────────────────────────────────────────────────────┘
```

The 418 TFs are the King mmc4 seed (`TF? != NA`); before the DE scan they are unioned with the unified master TF catalog (King ∪ Perez MOESM5), expanding the seed to **14,682 IDs** (`pipeline.py:189-227`). Only the 418-TF King subset is forced into the HVG panel (§1.9 item 9).
