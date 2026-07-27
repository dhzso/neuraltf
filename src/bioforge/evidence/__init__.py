"""BioForge Evidence Integration Framework (Layer 8B).

Cross-atlas transcription-factor prioritization for the NeuralTF project,
designed as a reusable framework (see ADR-0001 / ADR-0002).

Public API
----------
- :class:`EvidenceRecord`, :class:`ConfidenceTier`, :class:`EvidenceSource`
  — schema for per-TF evidence.
- :class:`BridgeTable`, :func:`load_bridge`, :func:`build_bridge_from_names`
  — gene identifier bridging (dd_Smed_v4 ↔ dd_Smed_v6).
- :class:`AtlasHarmonizer` — canonical tissue-label harmonization.
- :class:`EvidenceScorer`, :func:`rank_candidates` — multi-criterion scoring.
- :func:`assign_tiers`, :class:`ConfidencePolicy` — high/medium/low tiering.
- :func:`annotate_function` — minimal ontology mapping (stub).
- :mod:`bioforge.evidence.readers` — King/Fincher/Plass readers.
"""
from bioforge.evidence.confidence import ConfidencePolicy, assign_tiers
from bioforge.evidence.gene_mapping import BridgeTable, build_bridge_from_names, load_bridge
from bioforge.evidence.harmonization import CANONICAL_TISSUES, AtlasHarmonizer
from bioforge.evidence.ontology import annotate_function
from bioforge.evidence.scoring import DEFAULT_WEIGHTS, EvidenceScorer, rank_candidates
from bioforge.evidence.schema import ConfidenceTier, EvidenceRecord, EvidenceSource

__all__ = [
    "EvidenceRecord",
    "EvidenceSource",
    "ConfidenceTier",
    "BridgeTable",
    "load_bridge",
    "build_bridge_from_names",
    "AtlasHarmonizer",
    "CANONICAL_TISSUES",
    "EvidenceScorer",
    "DEFAULT_WEIGHTS",
    "rank_candidates",
    "ConfidencePolicy",
    "assign_tiers",
    "annotate_function",
]
