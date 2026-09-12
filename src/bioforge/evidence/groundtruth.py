"""Ground-truth labels for the King 2024 RNAi screen (mmc5).

PROBLEM THIS MODULE SOLVES (2026-09-11 review):
    King 2024 mmc5 is titled "All Transcription Factors Inhibited" — it
    lists every TF for which RNAi was PERFORMED and markers assayed, not
    only TFs with observed phenotypes. The table's own caption encodes
    phenotype status by font colour ("Transcription factors in red and
    marker genes in green have phenotypes"), but the distributed xlsx
    copy contains NO coloured runs (verified cell-by-cell: every cell is
    theme-1 black; sharedStrings carry no per-run colours). The
    pipeline's original `tested` label therefore meant "was in the
    screening list" — NOT "RNAi-validated".

    Reading the paper's main text + figure captions (Figs 3J, 4E, S4,
    S7A–H, S8A–G; every "RNAi showing loss of … following <TF> RNAi"
    statement) yields the definitive set of TFs with FISH-confirmed
    cell-type-loss phenotypes. Every listed phenotype is tied to a v6
    gene ID via the paper's own figure labels (e.g. "PHOX2A (dd29211)")
    or via the manual name->v6 map used by build_king_atlas.py
    (e.g. "Tbx2/3b" -> dd_Smed_v6_6470_0_1).

LABEL VOCABULARY (enforced here, consumed everywhere):
    screened             RNAi performed in King 2024 (in mmc5). Phenotype
                         status NOT implied.
    phenotype_confirmed  FISH-confirmed loss-of-cell-type phenotype shown
                         in the paper (Fig 3J/4E, S4, S7, S8). A strict
                         subset of `screened`.
    not_screened         No RNAi record in the King screen.

The `ProofStatus.TESTED` enum value keeps its string form ("tested") for
on-disk compatibility, but its DOCUMENTED meaning is now "screened".
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
#     Tbx2/3b  GLIPR1 (dd_210)+ loss   (4E, S8G)   -> dd_Smed_v6_6470_0_1
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
# dd17143 (TBX2) was screened but shows NO published phenotype; it is NOT
# in this set.
#
# COVERAGE NOTE: dd_Smed_v6_48508_0_1 (dd_1936+ neuron-loss phenotype, S8B)
# is phenotype-confirmed but is NOT part of the pipeline's candidate
# universe (its King-atlas rows are non-neural / below the neural gate and
# it was never seeded by another route) — an honest false negative of the
# discovery funnel, documented rather than force-seeded. Its RNAi row
# (dd48508, markers dd83746/dd1936) exists in mmc5, so a future candidate
# universe that seeds all mmc5 targets would include it.
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
    "dd_Smed_v6_6470_0_1",    # Tbx2/3b; GLIPR1+ loss
    "dd_Smed_v6_11500_0_1",   # IRX2; th+ (dopaminergic) loss
    "dd_Smed_v6_22163_0_1",   # UNCX; sert+ (serotonergic) loss
    "dd_Smed_v6_20282_0_1",   # dd_20282; dd_1248+ loss
    "dd_Smed_v6_22331_0_1",   # dd_22331; dd_29413+ loss
    "dd_Smed_v6_48508_0_1",   # dd_48508; dd_1936+ loss
    "dd_Smed_v6_11693_0_1",   # Tbx2/3c; dd_29413+ loss
    "dd_Smed_v6_8104_0_1",    # SOX2 / soxB1-2; dd_29413+ loss
})

# Short dd#### tokens for matching tables that only carry the numeric form.
# Used ONLY for the v4 dialect (isoform-blind by construction); v6 IDs are
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


def is_phenotype_confirmed(gene_id: str | None) -> bool:
    """True iff ``gene_id`` has a FISH-confirmed phenotype in King 2024.

    Matching discipline (mirrors the pipeline's 2026-09-06 isoform rule):
      - v6 IDs must match a curated v6 ID EXACTLY — a numeric-sibling
        isoform (e.g. dd_Smed_v6_10911_0_2 vs the curated _0_1) does NOT
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
    # Token fallback ONLY for the v4 dialect (v6 siblings stay distinct).
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
