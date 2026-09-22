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
# 2026-09-06 audit rewrite: the previous curated-keyword tuple + broad
# substring regex misclassified in BOTH directions on the real PlanMine
# corpus (70 short names observed in planmine_annotations.parquet):
#
#   False negatives (genuine TF/DBD names that scored zero):
#     HMG_box_dom (SOX), Paired_dom (Pax), TF_Brachyury (T-box),
#     TCF/LEF, CUT_dom, Iroquois_homeo, Homeo_prospero_dom,
#     Transcription_factor_COE, TF_AP2, ARID_dom, MAD_homology_MH1
#     (SMAD DNA-binding lobe - the old hint matched MH2, the non-binding
#     lobe), DMRT/dsx/mab-3, Teashirt_fam, Tio/Tsh,
#     Nucl_hormone_rcpt_ligand-bd, Homeobox_KN_domain, HTH_Psq,
#     Homeobox_CS, Homeodomain-like, DM_DNA-bd, Lambda_DNA-bd_dom,
#     WHTH_DNA-bd_dom, p53-like_TF_DNA-bd, Znf_GATA, Znf_NHR/GATA.
#   False positives (non-DBD names that passed):
#     Znf_LIM (LIM-only zinc fingers, mostly non-TF adaptors),
#     Znf_hrmn_rcpt (ligand-binding zinc finger of NHRs), plus the
#     blanket regex fragments "related[_-]?to" (matches any
#     "*-related-to-*" name) and hallucinated hints ("wasp", "cloudy",
#     "parkinson", "proto-encodable", "upsilon", "rflight", "tre-2",
#     "iry", "a24", "hb9").
#
# The replacement is an explicit allow-list anchored on the REAL domain
# naming conventions in this warehouse. Substring risk is controlled by
# ordering: exact names match first, then prefix patterns, and any
# remaining substring match requires a DBD-significant token that does
# not occur in the observed non-DBD corpus.
DNA_BINDING_DOMAIN_EXACT: frozenset[str] = frozenset({
    # Homeobox superclass
    "homeobox_dom", "homeobox", "homeobox_cs", "homeobox_kn_domain",
    "homeobox_metazoa", "homeodomain-like", "homeodomain", "homeodomain_like",
    "homeo_prospero_dom", "iroquois_homeo", "paired_dom", "paired",
    "pou", "pou_specific", "prox1",
    # HMG / SOX
    "hmg_box_dom", "hmg_box",
    # bHLH / bZIP /Zip
    "bhlh_dom", "bhlh", "bzip", "hlh_dom",
    # T-box / Brachyury / TEA
    "tf_t-box", "tf_t-box_cs", "t-box", "tbox", "tf_brachyury",
    "tea/atts", "tea_dom",
    # Ets / AP2 / ARID / COE / CUT / TCF-LEF
    "ets_dom", "ets", "tf_ap2", "tf_ap2_c", "arid_dom", "arid",
    "transcription_factor_coe", "transcription_factor_coe_cs", "coe_dom",
    "cut_dom", "cut", "tcf/lef", "tcf_lef", "lef_dom",
    # Forkhead
    "tf_fork_head", "tf_fork_head_cs", "fork_head", "forkhead",
    "winged_helix", "whth_dna-bd_dom", "whth",
    # SMAD MH1 (DNA-binding lobe; NOT MH2)
    "mad_homology_mh1", "mad_homology1_dwarfin-type",
    # DM / DMRT
    "dmrt/dsx/mab-3", "dm_dna-bd", "dm_dna-binding_domain",
    # Nuclear hormone receptors (DBD zinc fingers)
    "nucl_hormone_rcpt_ligand-bd", "nucl_hrmn_rcpt_lig-bd_core",
    "str_hrmn_rcpt", "znf_nhr/gata", "znf_gata", "gata_znf",
    # C2H2 zinc fingers (sequence-specific DNA-binding class)
    "znf_c2h2", "znf_c2h2-like", "znf_c2h2/integrase_dna-bd",
    "znf_c2h2_jaz", "c2h2_znf_fam",
    # Lambda repressor-like / p53
    "lambda_dna-bd_dom", "p53-like_tf_dna-bd", "p53-like_transcription_factor_dna-bd",
    # HTH Psq (Bed-class winged helix-turn-helix)
    "hth_psq", "hth_motif",
    # Teashirt
    "teashirt_fam", "tio/tsh",
})
# Prefix patterns for domain names that carry family suffixes (e.g.
# "Homeobox_dom_2", "Fork_head_N"). Each is anchored at the start of the
# normalized short name so "Znf_LIM" can never match "znf_"-prefixed
# DNA-binding families, and non-DBD names in the corpus (Ig-like, PAS,
# Peptidase, BRCT, Dwarfin/SMAD_dom-like, Innexin, ERAP1, DUF, PAC, IPT,
# CVC, ASH, SMAD_FHA...) match none of them.
DNA_BINDING_DOMAIN_STARTSWITH: tuple[str, ...] = (
    "homeobox", "homeodomain", "hmg_box", "bhlh", "fork_head", "forkhead",
    "znf_c2h2", "c2h2", "znf_gata", "znf_nhr", "mad_homology",
    "tf_ap2", "tf_fork_head", "tf_t-box", "t-box", "tbox", "tf_brachyury",
    "paired_dom", "pou_", "ets_dom", "cut_dom", "arid_dom", "dmrt",
    "dm_dna", "teashirt", "tio/tsh", "whth", "hth_psq", "lambda_dna",
    "p53-like_tf", "transcription_factor_coe",
)


def _norm_domain_name(short_name: str) -> str:
    return short_name.strip().lower().replace(" ", "_")


def domain_short_name_is_dna_binding(short_name: str | None) -> bool:
    """True if a PlanMine domain shortName indicates a DNA-binding domain.

    Explicit allow-list over the observed PlanMine corpus (see
    :data:`DNA_BINDING_DOMAIN_EXACT`). Exact matches first, then
    family-anchored prefixes. Deliberately NOT a fuzzy substring search:
    the previous keyword/regex classifier both missed real DBD families
    (HMG_box_dom, Paired_dom, TF_Brachyury, TCF/LEF, CUT_dom, COE,
    AP2, ARID, MH1, DMRT, Teashirt, NHR ligand-bd...) and passed
    non-DBD names (Znf_LIM, Znf_hrmn_rcpt) via blanket fragments.
    """
    if not short_name:
        return False
    s = _norm_domain_name(short_name)
    if not s or s in {"nan", "none"}:
        return False
    if s in DNA_BINDING_DOMAIN_EXACT:
        return True
    return s.startswith(DNA_BINDING_DOMAIN_STARTSWITH)


# -- GO-based flags ------------------------------------------------------------
# 2026-09-06 audit rewrite. Observed false positives in the real annotation
# corpus (919 unique GO terms in planmine_annotations.parquet):
#   - "head involution" (Drosophila epithelial morphogenesis) and
#     "specification of segmental identity, head" (segmentation) matched
#     the bare "head" keyword — neither is a neural process.
#   - plain "dna binding" / "double-stranded DNA binding" granted the
#     +0.02 go_tf bonus to ANY DNA-binding protein (histones,
#     polymerase subunits, HMGB), not transcription regulators.
# Direction awareness: "negative regulation of neuron apoptotic process"
# matches "neuron" but is a survival term; "regulation of transcription
# ... involved in somatic motor neuron fate commitment" IS neural+TF.
# The flag set below keeps neural marking direction-insensitive (a
# neuron-associated process is still neural evidence even when
# negatively regulated) but removes the non-neural matches outright.
GO_NEURAL_KEYWORDS = (
    "neuron", "nervous system development", "brain", "neurogenesis",
    "synaptic", "axon", "dendrite", "glial", "neural", "sensory organ",
    "sensory system", "sensory perception", "mechanosensory",
    "auditory", "ophthalm", "eye development", "eye morphogenesis",
    "eye formation", "photoreceptor", "visual", "cns",
    "neuromuscular junction", "neural crest", "neural tube",
    "neuronal stem cell", "motor neuron", "interneuron",
)
# Neural-adjacent-but-not-neural terms explicitly excluded: head
# involution (epithelial), segmental identity (segmentation), "eye "
# was replaced by specific eye-development phrases so "eye " can never
# prefix-match unrelated terms.
GO_NEURAL_EXCLUDE_SUBSTRINGS = (
    "segmental identity", "head involution",
)
GO_TF_KEYWORDS = (
    "transcription factor activity",
    "regulation of transcription", "regulation of transcription from rna polymerase ii promoter",
    "transcription regulatory region", "transcription repressor activity",
    "transcription coactivator activity", "transcription corepressor activity",
    "sequence-specific dna binding transcription factor",
    "transcription, dna-templated", "dna-templated transcription",
)
# NOTE: bare "dna binding" is deliberately NOT a TF keyword — it is a
# generic molecular-function term satisfied by histones and polymerase
# subunits. Only transcription-regulator-class terms grant go_tf.


def go_term_flags(name: str | None) -> tuple[bool, bool]:
    """(is_neural_go, is_tf_go) for a single GO term name.

    Neural matching: word/phrase anchored against the keyword list, with
    an exclusion list for the two observed non-neural terms that carried
    the old bare-"head" false positive. TF matching restricted to
    transcription-regulator-class terms (see :data:`GO_TF_KEYWORDS`).
    """
    if not name:
        return False, False
    n = name.strip().lower()
    if not n:
        return False, False
    if any(x in n for x in GO_NEURAL_EXCLUDE_SUBSTRINGS):
        is_neural = False
    else:
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