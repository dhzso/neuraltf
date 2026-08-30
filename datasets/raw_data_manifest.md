# NeuralTF Raw Data Manifest
*Verified Pipeline Audit & Dataset Inventory*

> **Note:** All processed outputs are generated deterministically from the raw downloads listed below.
> Raw files must never be modified. See `datasets/MANIFEST.md` for SHA-256 checksums and download URLs.

---

## 1. Active Primary Raw Files (Consumed Directly by Pipeline)

| File Name | Source Atlas | Data Type | Format & Size | Target Pipeline Script | Processed Destination |
|-----------|-------------|-----------|---------------|------------------------|-----------------------|
| `GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz` | Fincher (GSE111764) | Expression (whole-animal scRNA-seq) | .txt.gz ~96 MB | `scripts/convert_fincher.py` | `datasets/processed/fincher_subsample.h5ad` |
| `GSE103633_RAW.tar` | Plass (GSE103633) | Expression (11 samples scRNA-seq) | .tar ~26 MB | `scripts/consolidate_plass.py` | `datasets/processed/plass_v6.h5ad` |
| `adata_scRNA_Annotated.h5ad` | Cui (OMIX003867-01) | Expression (55,014 cells, 8 timepoints) | .h5ad ~1 GB | `projects/NeuralTF/scripts/convert_cui.py` | `datasets/processed/cui_v6.h5ad` |
| `1-s2.0-S2211124724001712-mmc4.xlsx` | King (2024) | Master TF Catalog (Gene ID, TF family, FSTF flag) | .xlsx ~80 KB | `scripts/build_master_catalog.py`, `pipeline.py` | `projects/NeuralTF/data/master_tf_catalog.csv` |
| `1-s2.0-S2211124724001712-mmc5.xlsx` | King (2024) | FSTF RNAi screen phenotypes | .xlsx ~12 KB | `pipeline.py` (RNAi stream), `prioritize_neural_tfs.py` | `runs/pipeline_run/rank.csv` |
| `1-s2.0-S2211124724001712-mmc6.xlsx` | King (2024) | TF pair correlations (X1/G0) | .xlsx ~10 KB | `pipeline.py` (Correlation stream) | `runs/pipeline_run/rank.csv` |
| `1-s2.0-S2211124724001712-mmc7.xlsx` | King (2024) | Neural G0 subcluster log2FC table | .xlsx ~96 KB | `scripts/build_king_atlas.py` | `projects/NeuralTF/data/king_atlas.tsv` |
| `41467_2025_65712_MOESM5_ESM.xlsx` | Perez (2025) | TF classification + h1SMcG↔v6 orthology | .xlsx ~16 MB | `scripts/preprocess_perez.py`, `pipeline.py` | `projects/NeuralTF/data/perez_tf_summary.csv` |
| `41467_2025_65712_MOESM19_ESM.xlsx` | Perez (2025) | ANANSE regulatory influence scores (9 fates) | .xlsx ~500 KB | `pipeline.py` (Perez influence stream) | `runs/pipeline_run/rank.csv` |
| `41467_2025_65712_MOESM22_ESM.xlsx` | Perez (2025) | ANANSE TF-target regulatory network (13,746 edges) | .xlsx ~859 KB | `projects/NeuralTF/scripts/ananse_full_scan.py` | `projects/NeuralTF/results/ananse_network_full.csv` |
| `smed_20140614.mapping.rosettastone.2020.txt` | PLANOSPHERE | Cross-assembly gene ID bridge (SMED ↔ v4 ↔ v6) | .txt ~67 MB | `scripts/build_bridge.py`, `convert_cui.py` | `projects/NeuralTF/data/bridge.csv` |
| `go.obo` | Gene Ontology | GO term definitions & namespaces | .obo ~37 MB | `projects/NeuralTF/scripts/make_supp_go_figures.py` | Supplementary figures |

---

## 2. Supplementary & Author-Deposited Reference Files (Unused by Core Pipeline)

The following files exist in `datasets/raw/` from upstream repository downloads or supplementary bundles, but are **not** consumed by the primary neural TF prioritization pipeline:

| File Name | Location | Type | Notes |
|-----------|----------|------|-------|
| `GSE111764_BrainClustering...` | `raw/GSE111764_GEO_Fincher_atlas/` | Expression (.txt.gz ~15 MB) | Supplementary brain-specific subset; whole-animal is used |
| `GSE111764_SexualClustering...` | `raw/GSE111764_GEO_Fincher_atlas/` | Expression (.txt.gz ~21 MB) | Supplementary sexual-cell subset |
| `GSE111764_Saturation.txt.gz` | `raw/GSE111764_GEO_Fincher_atlas/` | QC / Saturation | Reference metadata |
| `GSE103633_family.soft.gz` | `raw/GSE103633_GEO_Plass_atlas/` | GEO Metadata | Study annotations |
| `GSE103633_dd_Smed_v6.pcf.contigs.fasta.bz2` | `raw/GSE103633_GEO_Plass_atlas/` | Genome Fasta (~13 MB) | Assembly contigs |
| `adata_Neoblast.h5ad` | `raw/OMIX003867_OMIX_Cui_atlas/OMIX003867-01/` | H5AD (~200 MB) | Cui neoblast subset; full annotated dataset is used |
| `plk1_cut5d.h5ad` | `raw/OMIX003867_OMIX_Cui_atlas/OMIX003867-01/` | H5AD (~221 MB) | Cui perturbation dataset |
| `OMIX003867-02/Seurat_Visium/` | `raw/OMIX003867_OMIX_Cui_atlas/` | Spatial Visium (Robj) | Spatial transcriptomics; outside scRNA pipeline scope |
| `OMIX003867-03/Scanpy_Visium/` | `raw/OMIX003867_OMIX_Cui_atlas/` | Spatial Visium (H5AD) | Spatial transcriptomics; outside scRNA pipeline scope |
| `1-s2.0-S2211124724001712-mmc1.pdf` | `raw/Supplementary_Data_ King_2024/` | PDF | Author supplementary methods |
| `1-s2.0-S2211124724001712-mmc2.xlsx` | `raw/Supplementary_Data_ King_2024/` | Excel (~37 KB) | Cell cluster metadata |
| `1-s2.0-S2211124724001712-mmc3.xlsx` | `raw/Supplementary_Data_ King_2024/` | Excel (~308 KB) | Marker gene lists |
| `1-s2.0-S2211124724001712-mmc8.pdf` | `raw/Supplementary_Data_ King_2024/` | PDF | Author supplementary figures |
| `41467_2025_65712_MOESM*.pdf` | `raw/Supplementary_Data_ Perez_2025/` | PDFs | Author supplementary figures/notes |
| `41467_2025_65712_MOESM12/13_ESM.xlsx` | `raw/Supplementary_Data_ Perez_2025/` | Excel (~3 MB) | Lineage marker tables |

---

## 3. Atlas Processing Status Summary

| Atlas | Primary Paper | Accession / DOI | Cell / Gene Dimensions | Pipeline Processing Status |
|-------|---------------|-----------------|------------------------|----------------------------|
| **Fincher 2018** | Science 2018 | [GSE111764](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111764) | 50,587 cells × 28,069 genes (v4) | ✅ Processed (`fincher_subsample.h5ad`) |
| **Plass 2018** | Science 2018 | [GSE103633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103633) | 21,612 cells × 27,614 genes (v6) | ✅ Processed (`plass_v6.h5ad`) |
| **Cui 2023** | Nat Commun 2023 | [OMIX](https://ngdc.cncb.ac.cn/omix/release/OMIX003867) | 55,014 cells × 19,198 genes (v6) | ✅ Processed (`cui_v6.h5ad`, 2.05 GB) |
| **King 2024** | Cell Reports 2024 | [Cell Reports](https://www.sciencedirect.com/science/article/pii/S2211124724001712) | G0/X1 FACS single-cell TF atlas | ✅ Processed (`king_atlas.tsv`, mmc4–mmc7) |
| **Perez 2025** | Nat Commun 2025 | [Nature Communications](https://www.nature.com/articles/s41467-025-65712-0#Sec94) | TF classification + ANANSE GRN & influence | ✅ Processed (`perez_tf_summary.csv`, `ananse_network_full.csv`) |

