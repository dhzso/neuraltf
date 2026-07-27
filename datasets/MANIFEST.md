# Datasets Manifest

This directory holds the raw data and reference materials for BioForge's
First Consumer project, **NeuralTF** (planarian transcription factor
prioritization). The raw and reference files are **not tracked in git**
(see `.gitignore`). Reproduce the layout below by downloading from the
original sources or by copying from the project archive.

## Layout

```
datasets/
  raw/
    GSE103633_GEO_Plass_atlas/      # Plass et al., Science 2018
    GSE111764_GEO_Fincher_atlas/    # Fincher et al., Science 2018
    Supplementary_Data_King_2024/   # King et al., Cell Reports 2024
  reference/
    #Papers/                        # PDFs and summary of source papers
    Deepanshu_Master_slide.pptx     # MS thesis master slide
```

## Source-accession summary

| Dataset | Organism | Genome build | Cell atlas role |
|---------|----------|--------------|-----------------|
| GSE103633 (Plass 2018) | *Schmidtea mediterranea* | dd_Smed_v6 | Whole-anatomy atlas |
| GSE111764 (Fincher 2018) | *S. mediterranea* | dd_Smed_v4 | Brain/principal/sexual clusterings |
| King 2024 supplementary | *S. mediterranea* | dd_Smed_v6 | TF atlas + RNAi phenotypes (S4) + cluster log2FC (S6) |

## Gene-identifier bridge problem

Fincher uses **dd_Smed_v4** gene IDs; Plass and King use **dd_Smed_v6**. The
8B Evidence Integration Framework requires an explicit v4→v6 bridge table to
unify TF candidates across all three atlases. Do not guess by numeric prefix.
