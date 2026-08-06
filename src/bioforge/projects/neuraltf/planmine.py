"""PlanMine (InterMine) REST client for planarian functional annotations.

PlanMine is the community InterMine warehouse for planarian genomes and
transcriptomes.  *S. mediterranea* gene IDs of the ``dd_Smed_v6_*`` form are
stored as ``Contig`` records (a ``Transcript`` subclass) which carry:

- ``goAnnotation``         → GO Biological Process / Molecular Function
- ``domainHits``           → protein domain hits (Pfam / InterPro / SMART …)
- ``blastHits``            → best cross-species BLAST hits (incl. Human)
- ``sequence``             → full transcript nucleotide sequence

The host serves the legacy Tomcat chain without a full CA path, so this
client relaxes certificate verification deliberately (mirrors typical
PlanMine API usage).  Use strictly read-only, public, non-identifying data.

This module is import-time lean (``requests`` only) and fully unit-testable:
the query builders and the TF-domain / GO classifiers are pure functions.
"""

from __future__ import annotations

import logging
import re
import time
from typing import Any, Iterable

import requests
import urllib3

LOG = logging.getLogger(__name__)

PLANMINE_BASE_URL = "https://planmine.mpibpc.mpg.de/planmine/service"

# PlanMine's legacy Tomcat deployment serves a certificate chain that does not
# include the issuing CA, so standard verification fails.  The resource is
# read-only public data; relax verification for this request only and never
# reuse the flag for other hosts.
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# -- DNA-binding domain families (Pfam/InterPro short names) ------------------
# Curated keyword list used to decide whether a PlanMine protein-domain hit
# implies a "clear TF domain" for an uncharacterised candidate.  Matches by
# case-insensitive substring on the domain short name.
DNA_BINDING_DOMAIN_HINTS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("homeobox", ("homeobox", "homeodomain", "_hox", " hox ", "aristaless", "pou")),
    ("bHLH", ("basic-helix-loop-helix", "/bhlh", "hlh domain", "helix-loop-helix")),
    ("zinc-finger", ("znf_c2h2", "zf-c2h2", "c2h2-type", "c2h2_type", "c2h2/", "krueppel")),
    ("C4-type zinc finger", ("znf_c4", "c4-type")),
    ("forkhead", ("fork_head", "forkhead", "winged-helix")),
    ("HMG box", ("hmg box", "_sox", "hbp1")),
    ("bZIP", ("bzip", "basic leucine", "leucine zipper")),
    ("Ets", ("ets", "e_twenty", "wasp")),
    ("T-box", ("t-box", "_tbox", "tre-2")),
    ("MADS", ("mads", "minichromosome", "serum-response")),
    ("POU", ("pou", "iry")),
    ("SMAD", ("sman", "mh2", "wgr")),
    ("STAT", ("stat_", "/stat", "src_homology_2")),
    ("RUNT", ("_runt", "runt")),
    ("SAND", ("cloudy", "parkinson")),
    ("GCM", ("gcm", "a24")),
    ("TEA/ATTS", ("proto-encodable", "upsilon")),
    ("RFX", ("rflight", "_rfx")),
    ("Mothers/Caudal", ("caudal", "homothorax", "hb9")),
    ("Sox / SRY", ("sry", "_sox", "sx_")),
    ("Head-Activating", ("head_activation", "hairless-")),
)
# The list above keeps legitimate families; keep a canon 2nd list for a few
# high-signal DBD keywords that map to Pfam short names ending in patterns not
# covered by the first tuple (e.g. *_domain), add them in normalized lookup.
DNA_BINDING_DOMAIN_PREFIXES = (
    "homeobox", "homeodomain", "pou_", "hox", "fox_", "zf-", "znf-",
    "hlh", "sox_", "txn", "tf_", "myb_", "ets",
)


def domain_short_name_is_dna_binding(short_name: str | None) -> bool:
    """True if a PlanMine domain shortName indicates a DNA-binding domain."""
    if not short_name:
        return False
    s = short_name.strip().lower()
    if not s or s in {"nan", "none"}:
        return False
    for _, hints in DNA_BINDING_DOMAIN_HINTS:
        for hint in hints:
            if hint in s:
                return True
    return bool(re.search(
        r"(dna[_-]?bd|dna.binding|sequence.specific.dna|zinc[_-]?finger|"
        r"homeobox|homeodomain|forkhead|fork_head|helix.loop.helix|"
        r"bhlh|hlh|basic_region|bzip|pou|fn_winged_helix|hth[_-]|pta1|"
        r"myb|zf-c2h2|znf|c2h2|t[_-]?box|tbox|whth|schwann|related[_-]?to|"
        r"hox|sox|wg_shank)", s))


def match_dna_binding_family(short_name: str | None) -> list[str]:
    """Return the families (keys of :data:`DNA_BINDING_DOMAIN_HINTS`) a hit
    matches, best first."""
    if not short_name:
        return []
    s = short_name.strip().lower()
    if not s or s in {"nan", "none"}:
        return []
    out: list[str] = []
    for family, hints in DNA_BINDING_DOMAIN_HINTS:
        if any(h in s for h in hints):
            out.append(family)
    return out


# -- GO-based flags ------------------------------------------------------------
GO_NEURAL_KEYWORDS = (
    "neuron", "nervous system development", "brain", "neurogenesis",
    "synaptic", "axon", "dendrite", "glial", "neural", "sensory",
    "auditory", "ophthalm", "eye ", "visual", "head", "cns",
)
GO_TF_KEYWORDS = (
    "transcription factor activity", "dna binding", "dna-binding",
    "regulation of transcription", "nucleic acid",
)


def go_term_flags(name: str | None) -> tuple[bool, bool]:
    """(is_neural_go, is_tf_go) for a single GO term name."""
    if not name:
        return False, False
    n = name.strip().lower()
    is_neural = any(k in n for k in GO_NEURAL_KEYWORDS)
    is_tf = any(k in n for k in GO_TF_KEYWORDS)
    return is_neural, is_tf


class PlanMineError(RuntimeError):
    """Raised when PlanMine does not answer a query after retries."""


class PlanMineClient:
    """Small, resilient InterMine query client.

    Parameters
    ----------
    base_url : str
        PlanMine service URL (``.../planmine/service``).
    retries : int
        Maximum attempts per query (default 3).
    backoff : float
        Base delay (seconds) for exponential backoff between attempts.
    rate_limit : float
        Minimum pause (seconds) between requests to the service.
    logger : logging.Logger | None
        Logger for request/retry messages.
    """

    def __init__(
        self,
        base_url: str = PLANMINE_BASE_URL,
        *,
        retries: int = 3,
        backoff: float = 2.0,
        rate_limit: float = 0.25,
        logger: logging.Logger | None = None,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.retries = max(1, retries)
        self.backoff = backoff
        self.rate_limit = max(0.0, rate_limit)
        self._log = logger or LOG
        self._session: requests.Session | None = None

    # -- session ---------------------------------------------------------------
    @property
    def session(self) -> requests.Session:
        if self._session is None:
            self._session = requests.Session()
            self._log.info("PlanMine client target: %s", self.base_url)
        return self._session

    # -- low-level request -----------------------------------------------------
    def _get(self, endpoint: str, params: dict[str, Any]) -> requests.Response:
        url = f"{self.base_url}/{endpoint}"
        last_exc: Exception | None = None
        for attempt in range(1, self.retries + 1):
            try:
                resp = self.session.get(url, params=params, timeout=90, verify=False)
            except requests.RequestException as exc:  # network / TLS
                last_exc = exc
                self._log.warning(
                    "GET %s attempt %d/%d failed: %s",
                    endpoint, attempt, self.retries, type(exc).__name__,
                )
                time.sleep(self.backoff * attempt)
                continue
            if resp.status_code == 200:
                return resp
            # non-200: log and retry (transient 502/503, or 400 for bad XML
            # which is a coding bug — still retried harmlessly)
            self._log.warning(
                "GET %s attempt %d/%d -> HTTP %d",
                endpoint, attempt, self.retries, resp.status_code,
            )
            last_exc = PlanMineError(
                f"HTTP {resp.status_code} from {endpoint}: {resp.text[:80]!r}"
            )
            time.sleep(self.backoff * attempt)
        raise PlanMineError(f"PlanMine '{endpoint}' failed after retries") from last_exc

    # -- query layer ------------------------------------------------------------
    def query(
        self,
        view: Iterable[str],
        constraints: Iterable[tuple[str, str, str, str]],  # (path, op, value)
    ) -> list[list[str]]:
        """Run an InterMine query and return rows (tab-separated).

        ``view`` strings are model paths like ``Contig.primaryIdentifier``.
        Constraints are ``(path, op, value)`` triples AND-ed together.
        """
        xml = (
            f'<query model="genomic" view="{" ".join(view)}">'
            + "".join(
                f'<constraint path="{p}" op="{op}" value="{v}"/>'
                for (p, op, v) in constraints
            )
            + "</query>"
        )
        resp = self._get("query/results", {"query": xml})
        time.sleep(self.rate_limit)
        lines = resp.text.splitlines()
        return [l.split("\t") for l in lines if l.strip()]

    # -- contig (gene) annotation -----------------------------------------------
    def fetch_contig_annotations(self, contig_id: str) -> dict[str, Any]:
        """Fetch GO, protein domains, BLAST hits and sequence for one contig.

        Returns a dict with keys ``contig_id``, ``length``, ``sequence``,
        ``go_terms`` (list of dicts), ``domains`` (list of dicts),
        ``blast_hits`` (list of dicts).  Missing data = empty lists / None.
        """
        out: dict[str, Any] = {
            "contig_id": contig_id,
            "length": None,
            "sequence": None,
            "go_terms": [],
            "domains": [],
            "blast_hits": [],
        }
        con = (("Contig.primaryIdentifier", "=", contig_id),)

        # 1) GO terms
        for r in self.query(
            (
                "Contig.primaryIdentifier",
                "Contig.goAnnotation.ontologyTerm.name",
                "Contig.goAnnotation.ontologyTerm.namespace",
                "Contig.goAnnotation.ontologyTerm.identifier",
            ),
            con,
        ):
            if len(r) >= 4 and r[0] == contig_id:
                out["go_terms"].append(
                    {"name": r[1], "namespace": r[2], "identifier": r[3]}
                )

        # 2) protein domains (Pfam/InterPro short names)
        for r in self.query(
            (
                "Contig.primaryIdentifier",
                "Contig.domainHits.proteinDomain.shortName",
            ),
            con,
        ):
            if len(r) >= 2 and r[0] == contig_id:
                out["domains"].append({"short_name": r[1], "source": ""})

        # 3) cross-species BLAST hits (target, species, description)
        for r in self.query(
            (
                "Contig.primaryIdentifier",
                "Contig.blastHits.target",
                "Contig.blastHits.blastDomain.species",
                "Contig.blastHits.blastDomain.description",
            ),
            con,
        ):
            if len(r) >= 4 and r[0] == contig_id:
                out["blast_hits"].append(
                    {"target": r[1], "species": r[2], "description": r[3]}
                )

        # 4) length + full sequence
        for r in self.query(
            (
                "Contig.primaryIdentifier",
                "Contig.length",
                "Contig.sequence.residues",
            ),
            con,
        ):
            if len(r) >= 3 and r[0] == contig_id:
                try:
                    out["length"] = int(r[1])
                except (TypeError, ValueError):
                    out["length"] = None
                out["sequence"] = r[2]
        return out

    @staticmethod
    def human_ortholog(blast_hits: list[dict[str, Any]]) -> tuple[str, str] | None:
        """Best human BLAST hit as ``(species_label, description)`` or None."""
        for hit in blast_hits:
            sp = (hit.get("species") or "").strip().lower()
            if "homo sapiens" in sp or sp == "human":
                return hit["species"], hit["description"]
        return None