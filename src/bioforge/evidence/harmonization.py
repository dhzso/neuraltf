"""Atlas harmonization — align cell-type labels across atlases.

Fincher, Plass, and King each use their own cluster naming conventions.
The :class:`AtlasHarmonizer` maps each atlas's cluster labels onto a small
canonical tissue vocabulary so evidence integration can compare like with
like.

The canonical vocabularly is intentionally coarse (nine major tissue classes
matching the summary in ``datasets/reference``): ``epidermis``, ``muscle``,
``neuron``, ``intestine``, ``protonephridia``, ``cathepsin_positive``,
``parenchyma``, ``pharynx``, ``neoblast``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

from bioforge.core.logging import get_logger

logger = get_logger("evidence.harmonization")


CANONICAL_TISSUES: tuple[str, ...] = (
    "epidermis",
    "muscle",
    "neuron",
    "intestine",
    "protonephridia",
    "cathepsin_positive",
    "parenchyma",
    "pharynx",
    "neoblast",
)


@dataclass
class AtlasHarmonizer:
    """Map atlas-specific cluster labels to canonical tissue labels.

    Parameters
    ----------
    atlases
        Mapping of atlas name (e.g. ``"fincher"``) to a mapping of
        cluster-label → canonical-tissue. The framework supplies sensible
        defaults via :meth:`with_default_mappings`; callers can override.
    """

    atlases: dict[str, Mapping[str, str]] = field(default_factory=dict)

    def map(self, atlas: str, cluster_label: str) -> str:
        """Return the canonical tissue for a cluster label in `atlas`.

        Unknown labels fall back to ``"unknown"`` rather than raising, so
        downstream code can still score partial evidence.
        """
        try:
            mapping = self.atlases[atlas]
        except KeyError as exc:
            raise KeyError(f"unknown atlas '{atlas}'") from exc
        return mapping.get(cluster_label, "unknown")

    def add_atlas(self, atlas: str, mapping: Mapping[str, str]) -> None:
        self.atlases[atlas] = dict(mapping)
        logger.info("registered atlas '%s' with %d cluster labels",
                    atlas, len(mapping))

    @classmethod
    def with_default_mappings(cls) -> "AtlasHarmonizer":
        """Return a harmonizer seeded with the canonical tissue vocab.

        The default mappings are deliberately NOT populated with every
        cluster from every atlas — that is dataset-specific work. They
        simply pre-register the three atlas names with empty mappings
        so callers can extend them.
        """
        return cls(atlases={"fincher": {}, "plass": {}, "king": {}})
