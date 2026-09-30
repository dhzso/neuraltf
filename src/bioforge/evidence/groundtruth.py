"""Ground-truth labels for the King 2024 RNAi screen (mmc5).

King 2024 mmc5 is titled "All Transcription Factors Inhibited": it lists
every TF for which RNAi was performed and markers assayed, not only TFs with
observed phenotypes. Phenotype status is encoded by font colour in the
table's caption, but the distributed xlsx copy has no coloured runs, so the
FISH-confirmed subset here is read from the paper's own figures and captions
(Figs 3J, 4E, S4, S7A–H, S8A–G).

Label vocabulary (enforced here, consumed everywhere):
    screened             RNAi performed in King 2024 (in mmc5). Phenotype
                         status not implied.
    phenotype_confirmed  FISH-confirmed loss-of-cell-type phenotype shown
                         in the paper (Fig 3J/4E, S4, S7, S8). A strict
                         subset of `screened`.
    not_screened         No RNAi record in the King screen.

The `ProofStatus.TESTED` enum value keeps its string form ("tested") for
on-disk compatibility, but its documented meaning is "screened".
`phenotype_confirmed` is carried as a separate boolean column/annotation
so no downstream consumer can conflate the two.
"""
from __future__ import annotations

import re
from pathlib import Path

# --------------------------------------------------------------------------
# Curated phenotype-confirmed set (v6 IDs).
#
# Sources (King et al. 2024, Cell Reports 43:113843):
#   Parenchymal (Fig 3J, S4A–G):
#     dd_10911  glipr1+ loss            (S4A/S4C; also dd_6953+ neural loss, S7E)
#     Post-2b  dd_829 loss             (S4B)   -> dd_Smed_v6_9061_0_1
#     ascl-2   glipr1+ / dd_3069+ loss (S4C, 4E, S7F)
#     RREB2    mag-1+ / glipr1+ loss   (S4D)   -> dd_Smed_v6_10103_0_1
#     fer3l-1  SSPO (dd_628)+ loss     (S4E)   -> dd_Smed_v6_8096_0_1
#     IRX1     ZAN6 (dd_238)+ loss     (S4F)   -> dd_Smed_v6_14656_0_1
#     GCM2     KCP (dd_91)+ loss       (S4G)   -> dd_Smed_v6_7752_0_1
#   Neural (Fig 4E, S7A–H, S8A–G):
#     GFI1B    dd_28465+/dd_29413+ loss (4E, S7F)
#     IRX6     npp-18+ loss             (4E, S8A)  -> dd_Smed_v6_17352_0_1
#     INSM2    spp-4+ loss             (S7A)       -> dd_Smed_v6_28888_0_1
#     PHOX2A   dd_8060+ loss           (4E, S8E)   -> dd_Smed_v6_29211_0_1
#     POU4F3   CALM2 (dd_23127)+ loss  (4E, S8F)   -> dd_Smed_v6_30562_0_1
#     Tbx2/3b  GLIPR1 (dd_210)+ loss   (4E, S8G)   -> dd_Smed_v6_17143_0_1
#     IRX2     th+ loss                (S7F)       -> dd_Smed_v6_11500_0_1
#     UNCX     sert+ loss              (S7B)       -> dd_Smed_v6_22163_0_1
#     dd_20282 dd_1248+ loss           (S7C)
#     dd_22331 dd_29413+ loss          (S7D)
#     dd_48508 dd_1936+ loss           (S8B)
#     Tbx2/3c  dd_29413+ loss          (S8C)        -> dd_Smed_v6_11693_0_1
#     SOX2     dd_29413+ loss          (S8D)        -> dd_Smed_v6_8104_0_1
#     (soxB1-2, dd_8104, is additionally cited as the published regulator
#      of the CALM2+ population.)
#
# Tbx2/3b naming conflict: the paper's RNAi table (mmc5 row 69) lists
#     dd_17143 with marker dd_210 (GLIPR1) — the Tbx2/3b phenotype of Fig
#     4E/S8G — and King's mmc4 catalog names dd_17143 "T-box 2/3b protein"
#     (human TBX2 blast). dd_6470, which PlanMine labels "tbx2/3b", appears
#     in neither mmc5 nor mmc4, and its PlanMine descriptions are
#     aminopeptidase-related (alanyl aminopeptidase / CG14516) — an
#     internally inconsistent alias, most likely a Rosetta-Stone
#     many-to-many artifact. The ground truth follows the experiment:
#     Tbx2/3b = dd_Smed_v6_17143_0_1.
#
# COVERAGE NOTE: dd_Smed_v6_48508_0_1 (dd_1936+ neuron-loss phenotype, S8B)
# is phenotype-confirmed but is not part of the pipeline's candidate
# universe (its King-atlas rows are non-neural / below the neural gate and
# it was never seeded by another route) — an honest false negative of the
# discovery funnel, documented rather than force-seeded. Its RNAi row
# (dd48508, markers dd83746/dd1936) exists in mmc5, so a candidate universe
# that seeds all mmc5 targets would include it.
PHENOTYPE_CONFIRMED_V6: frozenset[str] = frozenset({
    "dd_Smed_v6_10911_0_1",   # zinc-finger TF; glipr1+ & dd_6953+ loss
    "dd_Smed_v6_9061_0_1",    # Post-2b; dd_829+ loss
    "dd_Smed_v6_14753_0_1",   # ascl-2; glipr1+ / dd_3069+ loss
    "dd_Smed_v6_10103_0_1",   # RREB2; mag-1+ / glipr1+ loss
    "dd_Smed_v6_8096_0_1",    # fer3l-1; SSPO+ loss
    "dd_Smed_v6_14656_0_1",   # IRX1; ZAN6+ loss
    "dd_Smed_v6_7752_0_1",    # GCM2; KCP+ loss
    "dd_Smed_v6_14824_0_1",   # GFI1B; dd_28465+/dd_29413+ loss
    "dd_Smed_v6_17352_0_1",   # IRX6; npp-18+ loss
    "dd_Smed_v6_28888_0_1",   # INSM2; spp-4+ loss
    "dd_Smed_v6_29211_0_1",   # PHOX2A; dd_8060+ loss
    "dd_Smed_v6_30562_0_1",   # POU4F3; CALM2+ loss
    "dd_Smed_v6_17143_0_1",   # Tbx2/3b; GLIPR1+ loss (the paper's actual
                               # RNAi target — see the naming note above)
    "dd_Smed_v6_11500_0_1",   # IRX2; th+ (dopaminergic) loss
    "dd_Smed_v6_22163_0_1",   # UNCX; sert+ (serotonergic) loss
    "dd_Smed_v6_20282_0_1",   # dd_20282; dd_1248+ loss
    "dd_Smed_v6_22331_0_1",   # dd_22331; dd_29413+ loss
    "dd_Smed_v6_48508_0_1",   # dd_48508; dd_1936+ loss
    "dd_Smed_v6_11693_0_1",   # Tbx2/3c; dd_29413+ loss
    "dd_Smed_v6_8104_0_1",    # SOX2 / soxB1-2; dd_29413+ loss
})

# Short dd#### tokens for matching tables that only carry the numeric form.
# Used only for the v4 dialect (isoform-blind by construction); v6 IDs are
# matched exactly so numeric-sibling isoforms never inherit the flag.
_RE_V6_NUM = re.compile(r"dd_Smed_v[46]_(\d+)_")
_RE_V4 = re.compile(r"dd_Smed_v4_(\d+)_")
PHENOTYPE_CONFIRMED_SHORT: frozenset[str] = frozenset(
    f"dd{_RE_V6_NUM.match(g).group(1)}" for g in PHENOTYPE_CONFIRMED_V6
)

# Gene-symbol aliases used in the paper for the same genes (Figure 4E/S7/S8
# label style). Used for name-form matching where dd#### is unavailable.
PHENOTYPE_CONFIRMED_NAMES: frozenset[str] = frozenset({
    "post-2b", "ascl-2", "rreb2", "fer3l-1", "irx1", "gcm2", "gfi1b",
    "irx6", "insm2", "phox2a", "pou4f3", "tbx2/3b", "irx2", "uncx",
    "tbx2/3c", "sox2", "soxb1-2",
})

# Curated symbol -> v6 map for the RNAi/phenotype genes that King's tables
# name symbolically (mmc5 stores 'GCM2', 'post2b', 'fer3l-1' with no dd####
# token). Every mapping is anchored to the paper's own figure labels (the
# PHENOTYPE_CONFIRMED_V6 curation above) or mmc4's GenBank descriptions,
# never guessed:
#   gcm2      Fig S4G KCP+ loss (dd7752)
#   post2b    Fig S4B dd_829+ loss; mmc5 'post2b (dd9061)' (dd9061)
#   fer3l-1   Fig S4E SSPO+ loss; mmc4 'fer3l-1 protein' (dd8096)
#   ascl-2    Fig S4C/4E glipr1+ / dd_3069+ loss (dd14753)
#   tbx2/3b   Fig 4E/S8G GLIPR1+ loss (dd17143, see the naming note above)
SYMBOL_TO_V6: dict[str, str] = {
    "gcm2": "dd_Smed_v6_7752_0_1",
    "post2b": "dd_Smed_v6_9061_0_1",
    "post-2b": "dd_Smed_v6_9061_0_1",
    "fer3l-1": "dd_Smed_v6_8096_0_1",
    "tbx2/3b": "dd_Smed_v6_17143_0_1",   # the paper's RNAi target, see above
    "ascl-2": "dd_Smed_v6_14753_0_1",
}


def is_phenotype_confirmed(gene_id: str | None) -> bool:
    """True iff ``gene_id`` has a FISH-confirmed phenotype in King 2024.

    Matching discipline (the same isoform rule the pipeline applies):
      - v6 IDs must match a curated v6 ID exactly — a numeric-sibling
        isoform (e.g. dd_Smed_v6_10911_0_2 vs the curated _0_1) does not
        inherit the flag; the paper's FISH panels name the gene, but the
        pipeline deliberately never aliases distinct isoforms.
      - v4-dialect IDs match via the bare dd#### token (the bridge maps
        v4 genes to v6 primary isoforms; v4 tables cannot distinguish
        isoforms anyway).
    Never guesses on any other form (names, junk, empty).
    """
    if not gene_id or not isinstance(gene_id, str):
        return False
    g = gene_id.strip()
    if g in PHENOTYPE_CONFIRMED_V6:
        return True
    # Token fallback only for the v4 dialect; v6 siblings stay distinct.
    if g.startswith("dd_Smed_v4_"):
        m = _RE_V4.match(g)
        if m and f"dd{m.group(1)}" in PHENOTYPE_CONFIRMED_SHORT:
            return True
    return False


def phenotype_status(rnai_score: float | None, gene_id: str | None) -> str:
    """Three-valued ground-truth status for a candidate.

    Returns one of:
      ``phenotype_confirmed``  in mmc5 AND FISH-confirmed phenotype
      ``screened_no_phenotype`` in mmc5, no phenotype shown in the paper
      ``not_screened``         no RNAi record
    """
    if rnai_score is not None and rnai_score > 0.0:
        return "phenotype_confirmed" if is_phenotype_confirmed(gene_id) \
            else "screened_no_phenotype"
    return "not_screened"


# --------------------------------------------------------------------------
# mmc5 ground-truth caveats shared by all consumers
# --------------------------------------------------------------------------
MMC5_GROUND_TRUTH_NOTE: str = (
    "King 2024 mmc5 lists ALL TFs inhibited by RNAi ('All Transcription "
    "Factors Inhibited'); phenotype status is encoded by font colour in the "
    "original table but the distributed xlsx copy is monochrome. 'tested' / "
    "rnai=1 therefore means SCREENED, not validated. "
    f"{len(PHENOTYPE_CONFIRMED_V6)} genes carry FISH-confirmed loss-of-"
    "cell-type phenotypes (Fig 3J/4E, S4, S7, S8 of the paper); see "
    "bioforge.evidence.groundtruth.PHENOTYPE_CONFIRMED_V6."
)

# Path helper so annotate scripts can locate this module's data if needed
MODULE_DIR = Path(__file__).resolve().parent
