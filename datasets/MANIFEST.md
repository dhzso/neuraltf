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
  references/
```

## Source-accession summary

| Dataset | Organism | Genome build | Cell atlas role |
|---------|----------|--------------|-----------------|
| GSE103633 (Plass 2018) | *Schmidtea mediterranea* | dd_Smed_v6 | Whole-anatomy atlas |
| GSE111764 (Fincher 2018) | *S. mediterranea* | dd_Smed_v4 | Brain/principal/sexual clusterings |
| King 2024 supplementary | *S. mediterranea* | dd_Smed_v6 | TF atlas + RNAi phenotypes (S4) + cluster log2FC (S6) |

## Source accessions

| Dataset | GEO | SRA | NCBI BioProject |
|---------|-----|-----|-----------------|
| Plass 2018 | [GSE103633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103633) | SRP117156 | PRJNA403817 |
| Fincher 2018 | [GSE111764](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111764) | SRP135258 | PRJNA438083 |
| King 2024 | Cell Reports [supplementary](https://www.sciencedirect.com/science/article/pii/S2211124724001712) (mmc4–mmc7) | — | — |

## Rosetta Stone gene-ID bridge

Fincher uses **dd_Smed_v4** gene IDs; Plass and King use **dd_Smed_v6**. The
v4↔v6 bridge table is built from the PLANOSPHERE Rosetta Stone mapping table
(2020), downloaded from:

<https://planosphere.stowers.org/pub/analysis/rosetta/smed_20140614.mapping.rosettastone.2020/smed_20140614.mapping.rosettastone.2020.txt>

Place it at `datasets/raw/smed_20140614.mapping.rosettastone.2020.txt`, then
run `python scripts/build_bridge.py`.

## Gene-identifier bridge problem

Fincher uses **dd_Smed_v4** gene IDs; Plass and King use **dd_Smed_v6**. The
8B Evidence Integration Framework requires an explicit v4→v6 bridge table to
unify TF candidates across all three atlases. Do not guess by numeric prefix.
