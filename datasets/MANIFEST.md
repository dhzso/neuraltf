# Datasets Manifest

Atlases processed artifacts are **generated
locally** from the raw downloads listed below — with a single command:

```bash
python scripts/generate_all.py            # everything, incl. PlanMine (network)
python scripts/generate_all.py --skip-planmine   # offline use
```

The generated outputs live in `datasets/processed/` (fincher/plass/cui h5ads,
PlanMine parquet + fasta), `projects/NeuralTF/data/` (bridge.csv,
king_atlas.tsv, master_tf_catalog.csv), `projects/NeuralTF/results/`
(Dirichlet CSVs, ANANSE network, supplementary tables S1–S4),
`projects/NeuralTF/figures/` (33 publication figures) and
`projects/NeuralTF/runs/pipeline_run/` (rank.csv, rank_neural.csv,
evidence_cards.md, 6 checkpoint parquets).

## Layout (download these to `datasets/raw/`): Keep the folder names same as below

```
datasets/
   raw/
    GSE103633_GEO_Plass_atlas/      # Plass et al., Science 2018 — GSE103633_RAW.tar, ~2 GB
    GSE111764_GEO_Fincher_atlas/    # Fincher et al., Science 2018 — DGE .txt.gz
    Supplementary_Data_King_2024/   # King et al., Cell Reports 2024 — mmc2–mmc7 xlsx
    OMIX003867_OMIX_Cui_atlas/      # Cui et al., 2023 — single‑cell h5ad (~55 GB) (download: [OMIX](https://ngdc.cncb.ac.cn/omix/release/OMIX003867))
      OMIX003867-01/
        singlecell_h5ad/
          adata_scRNA_Annotated.h5ad   # Cui annotated h5ad used by the pipeline
    smed_20140614.mapping.rosettastone.2020.txt    # Rosetta Stone gene-ID bridge, 67 MB
    go.obo                          # Gene Ontology, current release, ~40 MB
    Supplementary_Data_Perez_2025/   # Perez et al., 2025 — MOESM5 (TF classification), MOESM19 (ANANSE influence), MOESM22 (ANANSE GRN) (download: [Nature Communications](https://www.nature.com/articles/s41467-025-65712-0#Sec94))
  processed/                        # GENERATED: fincher_subsample.h5ad, plass_v6.h5ad, cui_v6.h5ad,
                                     # planmine_annotations.parquet, planmine_transcripts.fasta
  references/
```

## Source-accession summary

| Dataset | Organism | Genome build | Cell atlas role |
|---------|----------|--------------|-----------------|
| GSE103633 (Plass 2018) | *Schmidtea mediterranea* | dd_Smed_v6 | Whole-anatomy atlas |
| GSE111764 (Fincher 2018) | *S. mediterranea* | dd_Smed_v4 | Brain/principal/sexual clusterings |
| King 2024 supplementary | *S. mediterranea* | dd_Smed_v6 | TF atlas + RNAi phenotypes (S4) + cluster log2FC (S6) |
| Cui 2023 atlas | *S. mediterranea* | dd_Smed_v6 | Regeneration time‑course (8 time points, ~55 K cells |
| Perez 2025 supplementary | *S. mediterranea* | h1SMcG | TF classification (MOESM5), ANANSE influence (MOESM19), & ANANSE GRN validation (MOESM22) |
| Gene Ontology (go.obo) | model-organism-agnostic | dd_Smed_v6 | Canonical GO term names/namespaces for PlanMine annotation figures |

## Source accessions

| Dataset | GEO | SRA | NCBI BioProject / Accession |
|---------|-----|-----|-----------------|
| Plass 2018 | [GSE103633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103633) | SRP117156 | PRJNA403817 |
| Fincher 2018 | [GSE111764](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111764) | SRP135258 | PRJNA438083 |
| Cui 2023 | — | — | [OMIX](https://ngdc.cncb.ac.cn/omix/release/OMIX003867) |
| King 2024 | — | — | Cell Reports [supplementary](https://www.sciencedirect.com/science/article/pii/S2211124724001712) |
| Perez 2025 | — | — | Nature Communications [supplementary](https://www.nature.com/articles/s41467-025-65712-0#Sec94) |

## Data Integrity Checksums (SHA256)

**Verify after download:** `python scripts/verify_data.py`

| File | SHA256 | Size | Source |
|------|--------|------|--------|
| `datasets/raw/GSE103633_GEO_Plass_atlas/GSE103633_RAW.tar` | `FAC5E2B45BBF2E4CDFD0E1CE6241BAD3936C3EC667835869EA3DB6A4BD0B0DE2` | ~2 GB | GEO GSE103633 |
| `datasets/raw/GSE111764_GEO_Fincher_atlas/GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz` | `EF945898E17463D60BFED1C359A926383DD070A9098258599988640EBF0C29BC` | ~500 MB | GEO GSE111764 |
| `datasets/raw/Supplementary_Data_King_2024/1-s2.0-S2211124724001712-mmc4.xlsx` | `FD0F77946DCD6F35DCCC363761AB0857A4F60F94EBA4BA0282E289B6A5326A54` | ~5 MB | Cell Reports supplementary |
| `datasets/raw/Supplementary_Data_King_2024/1-s2.0-S2211124724001712-mmc5.xlsx` | `18A53AC54EFE6A244E875078D6863E67AC2B0E65D38DDA294615B05C3AB3E9D5` | ~1 MB | Cell Reports supplementary |
| `datasets/raw/Supplementary_Data_King_2024/1-s2.0-S2211124724001712-mmc6.xlsx` | `38646A0962CB577D719515CFA40A5B4E1569C5D3ADACA63DC51A37EDFFB2C6EF` | ~1 MB | Cell Reports supplementary |
| `datasets/raw/Supplementary_Data_King_2024/1-s2.0-S2211124724001712-mmc7.xlsx` | `49005C8E57810AE50DF95A6482036650CD75BF5206D9BBE49126364B9DE06B6F` | ~10 MB | Cell Reports supplementary |
| `datasets/raw/smed_20140614.mapping.rosettastone.2020.txt` | `D3E912B6798BD476B93C37E8B73F6801A1660C532D3EDD0E6BB8502ED14EC768` | 67 MB | PLANOSPHERE Rosetta Stone |
| `datasets/raw/go.obo` | `D3593751D885CA160B2AB7BAF6C7ECCD88CA3C4599F79436C674BAD661095FF0` | ~40 MB | Gene Ontology |
| `datasets/raw/Supplementary_Data_ Perez_2025/41467_2025_65712_MOESM5_ESM.xlsx` | `BB103659353D842E95487F85DDB89415DF978224A2088F6D203461055B7344FD` | ~1.6 MB | Perez 2025 supplementary |
| `datasets/raw/Supplementary_Data_ Perez_2025/41467_2025_65712_MOESM22_ESM.xlsx` | `9EDED529AC7E0D8AD0F5D6CD73B2102C59F3E829CC69DDB9200A530EADE5E458` | ~0.86 MB | Perez 2025 supplementary |
| `datasets/raw/OMIX003867_OMIX_Cui_atlas/OMIX003867-01/singlecell_h5ad/adata_scRNA_Annotated.h5ad` | `415B45A0A01C04DA616AA1909ADD605F05A13D514F1FAD72A2E622C5BD5DF414` | ~55 GB | Cui 2023 atlas annotated h5ad |

> **Note**: Verify checksums after download with `python scripts/verify_data.py` or `sha256sum <file>` (Linux/macOS) / `Get-FileHash -Algorithm SHA256 <file>` (Windows PowerShell).
> 
> **NOTE**: The pipeline integrates **five** atlases (Fincher, Plass, Cui, King, Perez) across **9 evidence streams**. Scoring weights: expression=0.200, all 8 others=0.100.

## Gene Ontology (go.obo)

Downloaded from the Gene Ontology current release — either the GO mirror at
<https://current.geneontology.org/ontology/go.obo> or the OBO Foundry URL
<http://purl.obolibrary.org/obo/go.obo> (plain text, ~40 MB). Place it at
`datasets/raw/go.obo`.

## Rosetta Stone gene-ID bridge

Fincher uses **dd_Smed_v4** gene IDs; Plass and King use **dd_Smed_v6**. The
v4↔v6 bridge table is built from the PLANOSPHERE Rosetta Stone mapping table
(2020), downloaded from:

<https://planosphere.stowers.org/pub/analysis/rosetta/smed_20140614.mapping.rosettastone.2020/smed_20140614.mapping.rosettastone.2020.txt>

Place the 67 MB table at
`datasets/raw/smed_20140614.mapping.rosettastone.2020.txt`, then run
`python scripts/build_bridge.py` (or `generate_all.py`) which writes
`projects/NeuralTF/data/bridge.csv`.

## Gene-identifier bridge problem

Fincher uses **dd_Smed_v4** gene IDs; Plass, Cui, and King use **dd_Smed_v6**. The Perez atlas uses the **h1SMcG** genome build whereafter pipeline maps between h1SMcG ↔ v6 via reciprocal best hit (RBH) orthology from MOESM5 [supplementary table from Perez et al 2025]. The
evidence integration framework requires an explicit **v6 bridge** table to unify TF candidates across all atlases.
