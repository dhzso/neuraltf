# Figure Guide — NeuralTF Publication Figures

## Visual Language

### Color Palette (colorblind-safe, Okabe-Ito inspired)

| Element | Color | Hex |
|---|---|---|
| Track A (RNAi-validated) | Strong blue | `#0072B2` |
| Track B (novel candidates) | Amber | `#E69F00` |
| Fixed-weight method | Near-black | `#333333` |
| Centered Dirichlet | Sky blue | `#56B4E9` |
| Uniform Dirichlet | Mauve | `#CC79A7` |
| Neural-filtered population | Mid-gray | `#999999` |
| Full 249 universe | Light-gray | `#DDDDDD` |
| Emphasis/highlight | Vermillion | `#D55E00` |

### Evidence Stream Colors

| Stream | Color |
|---|---|
| Expression | `#0072B2` |
| Specificity | `#E69F00` |
| Reproducibility | `#009E73` |
| RNAi | `#D55E00` |
| Correlation | `#CC79A7` |
| Neural enriched | `#56B4E9` |
| Neural specificity | `#F0E442` |

### Typography
- Font: Arial (fallback: Helvetica, DejaVu Sans)
- Base size: 7pt
- Panel labels: 9pt bold
- Axes titles: 8pt bold
- Axis labels: 7pt
- Tick labels: 6pt

### Dimensions
- Single column: 89mm (3.5 in)
- Double column: 180mm (7.08 in)
- DPI: 300 (PNG)
- Output: PNG + PDF

---

## Main Figures

### Figure 1 — Candidate prioritization landscape

**Scientific question:** How is the candidate space progressively narrowed from all TFs to the final prioritized candidates?

**Panels:**
- **A** Filtering flow: horizontal bar chart showing population at each stage (249 → expression-supported → neural-filtered → Track A → Track B → Top 10)
- **B** Score ECDF: cumulative distribution of integrated scores for all 249, with neural-filtered and Top 10 highlighted
- **C** Ranked landscape: scatter of rank vs integrated score with final candidates annotated

**Data sources:** `rank.csv` (249), `rank_neural.csv` (99), `top10_neural_tfs_prioritized.csv`

**Key statistics:** Population counts at each filtering stage, cumulative score distributions

**Intended conclusion:** The filtering pipeline progressively narrows from 249 to 10 candidates, with the final Top 10 occupying the top of the score distribution.

---

### Figure 2 — Evidence architecture of prioritized candidates

**Scientific question:** What evidence patterns distinguish high-priority candidates?

**Panels:**
- **A** Evidence heatmap: Top 10 × 7 streams with metadata tracks (proof status, neural specificity), cell values showing raw evidence strength
- **B** Score decomposition: stacked horizontal bars showing weighted contribution of each stream to the integrated score

**Data sources:** `rank_neural.csv`, `top10_neural_tfs_prioritized.csv`

**Key statistics:** Raw evidence values, weighted contributions, integrated scores

**Intended conclusion:** Each highly-ranked candidate has a distinct evidence profile, with different streams contributing to their prioritization.

---

### Figure 3 — Ranking robustness under weight uncertainty

**Scientific question:** How robust are the candidate rankings to assumptions about evidence-stream weights?

**Panels:**
- **A** Bump chart: parallel-rank visualization showing how Top 10 candidates rank under fixed, centered Dirichlet, and uniform Dirichlet weighting
- **B** Rank uncertainty: uncertainty intervals from 1000 weight sensitivity draws, with P(Top10) annotations
- **C** Robustness classification: candidates classified as stable (≥80% in Top 10), moderate (50-80%), or sensitive (<50%)

**Data sources:** `top10_neural_tfs_prioritized.csv`, `dirichlet_top10_prioritized.csv`, `dirichlet_uniform_top10.csv`, `weight_sensitivity_draws.csv`, `weight_sensitivity_top10_challengers.csv`

**Key statistics:** Rank positions across methods, rank distributions from simulations, P(Top10)

**Intended conclusion:** Core candidates (dd14115, dd19890, dd26877, dd38342, dd13343) are highly robust to weighting assumptions; some candidates are sensitive to the choice of prior.

---

### Figure 4 — Sensitivity of evidence streams

**Scientific question:** Which evidence streams are most important for candidate prioritization?

**Panels:**
- **A** Global impact: lollipop plot of median |rank change| when each stream is removed, with annotations of Top 10 displacement count
- **B** Candidate sensitivity: heatmap showing rank change per candidate × stream removed

**Data sources:** `rank_neural.csv`, `top10_neural_tfs_prioritized.csv`

**Key statistics:** Median rank change per stream removal, per-candidate rank shifts, number of Top 10 changes

**Intended conclusion:** RNAi and expression streams have the largest global impact; individual candidates show variable sensitivity to specific streams.

---

### Figure 5 — Neural filtering vs full candidate universe

**Scientific question:** Do the same candidates remain prioritized in the full 249-candidate universe?

**Panels:**
- **A** Rank-rank comparison: scatter of 99-neural rank vs 249-wide rank (uniform Dirichlet), with Spearman correlation
- **B** Score ECDF: cumulative score distributions for 99 neural vs 249 full, with KS test

**Data sources:** `rank.csv`, `rank_neural.csv`, `top10_neural_tfs_prioritized.csv`, `dirichlet_uniform_full_rank.csv`, `dirichlet_uniform_all249_full_rank.csv`

**Key statistics:** Spearman rank correlation, KS statistic and p-value

**Intended conclusion:** The neural-filtered top candidates are robustly recovered in the full 249-candidate universe (Spearman ρ > 0.9).

---

### Figure 6 — Prioritized candidate atlas

**Scientific question:** What is the complete evidence profile for each prioritized candidate?

**Panels:**
- **A** Candidate cards: compact visual encoding of gene identity, track, composite score, proof status, evidence stream dots, domain annotations, and human orthologs
- **B** Cross-method comparison: lollipop plot comparing composite scores across fixed, centered, and uniform methods

**Data sources:** `rank_neural.csv`, `top10_neural_tfs_prioritized.csv`, `dirichlet_top10_prioritized.csv`, `dirichlet_uniform_top10.csv`

**Key statistics:** Composite scores across three methods, proof status, domain annotations

**Intended conclusion:** The Top 10 candidates are supported by multiple lines of evidence and maintain high rankings across all weighting methods.

---

## Supplementary Figures

| Figure | Question | Data |
|---|---|---|
| S1: GO gene-term map | Which GO terms annotate neural TF candidates? | PlanMine annotations + go.obo |
| S2: GO top-10 dot matrix | Which GO terms are enriched in top candidates? | Same |
| S3: GO top-term statistics | Distribution of GO term annotations | Same |
| S4: GO neural focus | Neural-specific GO terms only | Same |

---

## Regeneration

```bash
# All main figures
python projects/NeuralTF/scripts/generate_publication_figures.py

# Specific figure
python projects/NeuralTF/scripts/generate_publication_figures.py --figure 1 3

# Individual figure
python projects/NeuralTF/scripts/figures/Fig1_candidate_landscape.py
```
