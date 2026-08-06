# Candidate Summary Report — NeuralTF Prioritization

Inputs: `rank_neural.csv` (96 neural candidates), PlanMine annotations (`datasets/processed/planmine_annotations.parquet`), the v6→v4 identifier bridge, King 2024 supplementary tables (mmc4 TF catalog, mmc5 FSTF RNAi screen), G0 atlas (`king_atlas.tsv`).

## Method

Two independent tracks:

- **Track A** — `proof_status == known_rnai_validated`: RNAi-validated benchmark TFs from the King 2024 FSTF screen. Top 5 by composite score.
- **Track B** — `proof_status == novel_candidate`: no published RNAi data; filtered to candidates with a clear DNA-binding TF domain (PlanMine protein-domain hits or mmc4 TF flag), then top 5 by composite score.

`composite_score = integrated_score + bonuses` (formula in `bioforge/projects/neuraltf/prioritize.py`): TF domain +0.05, neural GO +0.03, TF GO +0.02, human ortholog +0.02, RNAi-validated +0.02.

## Shortlist

| v6 id | gene_name | track | rank | composite | human ortholog |
|---|---|---|---|---|---|
| dd_Smed_v6_19890_0_1 | dd19890 | A | 1 | 0.937 | T-box transcription factor TBX20 isoform 1 [Homo sapiens] |
| dd_Smed_v6_2946_0_1 | dd2946 | A | 2 | 0.923 | zinc finger protein Aiolos isoform 3 [Homo sapiens] |
| dd_Smed_v6_26877_0_1 | dd26877 | A | 3 | 0.921 | protein atonal homolog 1 [Homo sapiens] |
| dd_Smed_v6_14115_0_1 | dd14115 | A | 4 | 0.918 | LIM homeobox transcription factor 1-alpha isoform X2 [Homo sapiens] |
| dd_Smed_v6_38342_0_1 | dd38342 | A | 5 | 0.905 | POU domain, class 3, transcription factor 4 [Homo sapiens] |
| dd_Smed_v6_31217_0_1 | dd31217 | B | 1 | 0.797 | helix-loop-helix protein 1 [Homo sapiens] |
| dd_Smed_v6_12170_0_1 | dd12170 | B | 2 | 0.736 | forkhead box protein J2 [Homo sapiens] |
| dd_Smed_v6_11930_0_1 | dd11930 | B | 3 | 0.731 |  |
| dd_Smed_v6_15253_0_1 | pitx | B | 4 | 0.712 | pituitary homeobox 3 [Homo sapiens] |
| dd_Smed_v6_14753_0_1 | ascl-2 | B | 5 | 0.709 | achaete-scute homolog 3 [Homo sapiens] |

## dd19890 (Track A, rank 1)

- v6 `dd_Smed_v6_19890_0_1` · v4 `dd_Smed_v4_19890_0_1` · `known_rnai_validated`
- composite `0.937` (pipeline integrated `0.797`, 6 evidence streams)
- DNA-binding domains (PlanMine): `TF_T-box, TF_T-box_CS, p53-like_TF_DNA-bd`
- GO terms: DNA binding; RNA polymerase II activating transcription factor binding; RNA polymerase II regulatory region sequence-specific DNA binding; RNA polymerase II transcription coactivator activity; RNA polymerase II transcription factor binding; aortic valve development; aortic valve morphogenesis; atrial septum morphogenesis; blood circulation; cardiac chamber formation; cardiac muscle tissue morphoge
- Human ortholog: T-box transcription factor TBX20 isoform 1 [Homo sapiens]
- RNAi note: RNAi-validated (King 2024, mmc5)
- Cross-stage dynamics: G0 progenitor max log2FC `5.59` · X1 n/a
- Wet-lab suggestion: design dsRNA against nt 300–800 of `dd_Smed_v6_19890_0_1` in `datasets/processed/planmine_transcripts.fasta`; FISH probe ≈ 800 nt antisense over the CDS region.

## dd2946 (Track A, rank 2)

- v6 `dd_Smed_v6_2946_0_1` · v4 `dd_Smed_v4_2946_0_1` · `known_rnai_validated`
- composite `0.923` (pipeline integrated `0.833`, 6 evidence streams)
- DNA-binding domains (PlanMine): `Znf_C2H2, Znf_C2H2-like, Znf_C2H2/integrase_DNA-bd`
- Human ortholog: zinc finger protein Aiolos isoform 3 [Homo sapiens]
- RNAi note: RNAi-validated (King 2024, mmc5)
- Cross-stage dynamics: G0 progenitor max log2FC `6.14` · X1 n/a
- Wet-lab suggestion: design dsRNA against nt 300–800 of `dd_Smed_v6_2946_0_1` in `datasets/processed/planmine_transcripts.fasta`; FISH probe ≈ 800 nt antisense over the CDS region.

## dd26877 (Track A, rank 3)

- v6 `dd_Smed_v6_26877_0_1` · v4 `dd_Smed_v4_26877_0_1` · `known_rnai_validated`
- composite `0.921` (pipeline integrated `0.781`, 6 evidence streams)
- DNA-binding domains (PlanMine): `bHLH_dom`
- GO terms: DNA binding; apoptotic process; auditory receptor cell differentiation; auditory receptor cell fate determination; auditory receptor cell fate specification; axon guidance; brain development; cell differentiation; central nervous system development; cerebral cortex development; chromatin DNA binding; inner ear development; inner ear morphogenesis; multicellular organismal development; nervous syst
- Human ortholog: protein atonal homolog 1 [Homo sapiens]
- RNAi note: RNAi-validated (King 2024, mmc5)
- Cross-stage dynamics: G0 progenitor max log2FC `5.00` · X1 n/a
- Wet-lab suggestion: design dsRNA against nt 300–800 of `dd_Smed_v6_26877_0_1` in `datasets/processed/planmine_transcripts.fasta`; FISH probe ≈ 800 nt antisense over the CDS region.

## dd14115 (Track A, rank 4)

- v6 `dd_Smed_v6_14115_0_1` · v4 `dd_Smed_v4_14115_0_1` · `known_rnai_validated`
- composite `0.918` (pipeline integrated `0.778`, 6 evidence streams)
- DNA-binding domains (PlanMine): `Homeobox_CS, Homeobox_dom, Homeodomain-like, Znf_LIM`
- GO terms: DNA binding; axon guidance; camera-type eye development; cell death; cell proliferation; central nervous system neuron development; central nervous system neuron differentiation; cerebellum development; cerebellum morphogenesis; collagen fibril organization; dentate gyrus development; dopaminergic neuron differentiation; embryonic limb morphogenesis; limb morphogenesis; metal ion binding; midbrain
- Human ortholog: LIM homeobox transcription factor 1-alpha isoform X2 [Homo sapiens]
- RNAi note: RNAi-validated (King 2024, mmc5)
- Cross-stage dynamics: G0 progenitor max log2FC `7.82` · X1 n/a
- Wet-lab suggestion: design dsRNA against nt 300–800 of `dd_Smed_v6_14115_0_1` in `datasets/processed/planmine_transcripts.fasta`; FISH probe ≈ 800 nt antisense over the CDS region.

## dd38342 (Track A, rank 5)

- v6 `dd_Smed_v6_38342_0_1` · v4 `dd_Smed_v4_38342_0_1` · `known_rnai_validated`
- composite `0.905` (pipeline integrated `0.765`, 6 evidence streams)
- DNA-binding domains (PlanMine): `Homeobox_CS, Homeobox_dom, Homeodomain-like, Lambda_DNA-bd_dom, POU, POU_specific`
- GO terms: AT DNA binding; DNA binding; cochlea morphogenesis; double-stranded DNA binding; forebrain neuron differentiation; inner ear development; intracellular part; negative regulation of mesenchymal cell apoptotic process; nucleus; regulation of transcription, DNA-templated; sensory perception of sound; sequence-specific DNA binding; transcription factor activity, sequence-specific DNA binding; transcri
- Human ortholog: POU domain, class 3, transcription factor 4 [Homo sapiens]
- RNAi note: RNAi-validated (King 2024, mmc5)
- Cross-stage dynamics: G0 progenitor max log2FC `4.40` · X1 n/a
- Wet-lab suggestion: design dsRNA against nt 300–800 of `dd_Smed_v6_38342_0_1` in `datasets/processed/planmine_transcripts.fasta`; FISH probe ≈ 800 nt antisense over the CDS region.

## dd31217 (Track B, rank 1)

- v6 `dd_Smed_v6_31217_0_1` · v4 `dd_Smed_v4_31217_0_1` · `novel_candidate`
- composite `0.797` (pipeline integrated `0.677`, 5 evidence streams)
- DNA-binding domains (PlanMine): `bHLH_dom`
- GO terms: DNA binding; cell differentiation; central nervous system development; multicellular organismal development; nucleus; protein dimerization activity; regulation of transcription, DNA-templated; transcription, DNA-templated
- Human ortholog: helix-loop-helix protein 1 [Homo sapiens]
- RNAi note: Not RNAi-tested in King 2024 mmc5; novel neural-fate candidate
- Cross-stage dynamics: G0 progenitor max log2FC `7.71` · X1 n/a
- Wet-lab suggestion: design dsRNA against nt 300–800 of `dd_Smed_v6_31217_0_1` in `datasets/processed/planmine_transcripts.fasta`; FISH probe ≈ 800 nt antisense over the CDS region.

## dd12170 (Track B, rank 2)

- v6 `dd_Smed_v6_12170_0_1` · v4 `dd_Smed_v4_12170_0_1` · `novel_candidate`
- composite `0.736` (pipeline integrated `0.646`, 5 evidence streams)
- DNA-binding domains (PlanMine): `TF_fork_head, TF_fork_head_CS, WHTH_DNA-bd_dom`
- GO terms: DNA binding; DNA binding, bending; double-stranded DNA binding; embryo development; nucleolus; nucleus; organ development; pattern specification process; positive regulation of transcription, DNA-templated; protein domain specific binding; regulation of sequence-specific DNA binding transcription factor activity; regulation of transcription, DNA-templated; sequence-specific DNA binding; tissue dev
- Human ortholog: forkhead box protein J2 [Homo sapiens]
- RNAi note: Not RNAi-tested in King 2024 mmc5; novel neural-fate candidate
- Cross-stage dynamics: G0 progenitor max log2FC `4.33` · X1 n/a
- Wet-lab suggestion: design dsRNA against nt 300–800 of `dd_Smed_v6_12170_0_1` in `datasets/processed/planmine_transcripts.fasta`; FISH probe ≈ 800 nt antisense over the CDS region.

## dd11930 (Track B, rank 3)

- v6 `dd_Smed_v6_11930_0_1` · v4 `dd_Smed_v4_11930_0_1` · `novel_candidate`
- composite `0.731` (pipeline integrated `0.661`, 5 evidence streams)
- DNA-binding domains (PlanMine): `Znf_C2H2, Znf_C2H2-like`
- RNAi note: Not RNAi-tested in King 2024 mmc5; novel neural-fate candidate
- Cross-stage dynamics: G0 progenitor max log2FC `7.11` · X1 n/a
- Wet-lab suggestion: design dsRNA against nt 300–800 of `dd_Smed_v6_11930_0_1` in `datasets/processed/planmine_transcripts.fasta`; FISH probe ≈ 800 nt antisense over the CDS region.

## pitx (Track B, rank 4)

- v6 `dd_Smed_v6_15253_0_1` · v4 `dd_Smed_v4_15253_0_1` · `novel_candidate`
- composite `0.712` (pipeline integrated `0.592`, 6 evidence streams)
- DNA-binding domains (PlanMine): `Homeobox_CS, Homeobox_dom, Homeodomain-like`
- GO terms: DNA binding; dendrite morphogenesis; dopaminergic neuron differentiation; lens development in camera-type eye; lens fiber cell differentiation; lens morphogenesis in camera-type eye; locomotory behavior; midbrain development; multicellular organismal development; neuron development; nucleus; organ morphogenesis; positive regulation of transcription, DNA-templated; protein binding; regulation of ge
- Human ortholog: pituitary homeobox 3 [Homo sapiens]
- RNAi note: Not RNAi-tested in King 2024 mmc5; novel neural-fate candidate
- Cross-stage dynamics: G0 progenitor max log2FC `6.26` · X1 n/a
- Wet-lab suggestion: design dsRNA against nt 300–800 of `dd_Smed_v6_15253_0_1` in `datasets/processed/planmine_transcripts.fasta`; FISH probe ≈ 800 nt antisense over the CDS region.

## ascl-2 (Track B, rank 5)

- v6 `dd_Smed_v6_14753_0_1` · v4 `dd_Smed_v4_14753_0_1` · `novel_candidate`
- composite `0.709` (pipeline integrated `0.619`, 5 evidence streams)
- DNA-binding domains (PlanMine): `bHLH_dom`
- GO terms: DNA binding; intracellular membrane-bounded organelle; nucleolus; nucleus; protein dimerization activity; regulation of transcription from RNA polymerase II promoter; transcription factor complex; transcription, DNA-templated
- Human ortholog: achaete-scute homolog 3 [Homo sapiens]
- RNAi note: Not RNAi-tested in King 2024 mmc5; novel neural-fate candidate
- Cross-stage dynamics: G0 progenitor max log2FC `6.66` · X1 n/a
- Wet-lab suggestion: design dsRNA against nt 300–800 of `dd_Smed_v6_14753_0_1` in `datasets/processed/planmine_transcripts.fasta`; FISH probe ≈ 800 nt antisense over the CDS region.

## Reproducibility

Deterministic pipeline (seed 42, pinned inputs, commit-pinned raw files). PlanMine snapshot dated on run; identifier mapping via the bridge table with explicit ambiguity flags (no numeric guessing).