"""Dual-track prioritization of NeuralTF candidates.

Combines the pipeline ranking (``rank_neural.csv``), PlanMine functional
annotations, the v6<->v4 identifier bridge, and the King 2024 supplementary
tables into a transparent, reproducible shortlist: top 5 RNAi-validated
benchmark TFs (Track A) and top 5 uncharacterized novel TFs (Track B).

Scoring
-------
``composite_score`` (0-1) is the pipeline ``integrated_score`` plus small,
fully documented bonuses:

+----------------------+------+------------------------------------------------+
| bonus                | wt   | conditions                                    |
+----------------------+------+------------------------------------------------+
| tf_domain            | 0.05 | DNA-binding domain (PlanMine) or TF flag in   |
|                      |      | the King mmc4 catalog                         |
+----------------------+------+------------------------------------------------+
| go_neural            | 0.03 | any GO term flagged neural                    |
+----------------------+------+------------------------------------------------+
| go_tf                | 0.02 | any GO term flagged transcription-related     |
+----------------------+------+------------------------------------------------+
| human_ortholog       | 0.02 | clear human BLAST ortholog (PlanMine or mmc4) |
+----------------------+------+------------------------------------------------+
| brain_phenotype (A)  | 0.02 | RNAi screen phenotype in King mmc5            |
+----------------------+------+------------------------------------------------+

Bonuses are small and additive; the pipeline score stays dominant.

All selector logic is pure and unit-testable (no I/O).
"""

from __future__ import annotations

import re
from typing import Any

import pandas as pd

from bioforge.projects.neuraltf.planmine import (
    domain_short_name_is_dna_binding,
    go_term_flags,
)

# bonus weights (see docstring)
BONUS_TF_DOMAIN = 0.05
BONUS_GO_NEURAL = 0.03
BONUS_GO_TF = 0.02
BONUS_HUMAN_ORTHOLOG = 0.02
BONUS_BRAIN_PHENOTYPE = 0.02

MAX_COMPOSITE = 1.0


# ---------------------------------------------------------------------------
# identifier mapping (v6 <-> v4)
# ---------------------------------------------------------------------------
def map_v6_to_v4(bridge: pd.DataFrame) -> pd.DataFrame:
    """One row per unique v6 id with its v4 id(s).

    When a v6 id maps to more than one distinct v4 id the mapping is flagged
    ``ambiguous`` and the v4 value is left blank - we never guess numeric
    prefix matches (the documented behaviour).
    """
    col_map: dict[str, str] = {}
    for c in bridge.columns:
        lc = str(c).strip().lower()
        if "v6" in lc:
            col_map.setdefault("v6", c)
        elif "v4" in lc:
            col_map.setdefault("v4", c)
    if "v6" not in col_map or "v4" not in col_map:
        raise ValueError(f"bridge table needs v6/v4 columns, got {list(bridge.columns)}")
    v6_col, v4_col = col_map["v6"], col_map["v4"]
    g = bridge[[v6_col, v4_col]].dropna(subset=[v6_col]).copy()
    g[v6_col] = g[v6_col].astype(str).str.strip()
    g[v4_col] = g[v4_col].astype(str).str.strip()
    g = g[g[v4_col].notna() & (g[v4_col] != "nan") & (g[v4_col] != "")]

    rows = []
    for v6, sub in g.groupby(v6_col, sort=False):
        v4_vals = sorted({v for v in sub[v4_col] if v})
        if len(v4_vals) == 1:
            rows.append({"v6_id": v6, "v4_id": v4_vals[0], "mapping_flag": "unique"})
        elif v4_vals:
            rows.append({"v6_id": v6, "v4_id": "", "mapping_flag": "ambiguous"})
        else:
            rows.append({"v6_id": v6, "v4_id": "", "mapping_flag": "unmapped"})
    out = pd.DataFrame(rows, columns=["v6_id", "v4_id", "mapping_flag"])
    return out if out.empty else out.sort_values("v6_id").reset_index(drop=True)


def attach_v4(rank: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    """Join v6->v4 mapping to the rank table, flagging ambiguities."""
    out = rank.merge(
        mapping, left_on="gene_id", right_on="v6_id", how="left"
    )
    out["gene_id_v4"] = out["v4_id"].fillna("")
    out["v4_mapping_flag"] = out["mapping_flag"].fillna("unmapped")
    return out.drop(columns=["v4_id", "mapping_flag", "v6_id"], errors="ignore")


# ---------------------------------------------------------------------------
# feature extraction
# ---------------------------------------------------------------------------
def extract_gene_symbol(describe: str) -> str:
    """Pull a gene symbol out of an NCBI protein title.

    e.g. ``"ALX homeobox protein 1 isoform X1 [Homo sapiens]"`` -> ``ALX``,
    ``"LIM/homeobox protein Lhx6 isoform 6 [Homo sapiens]"`` -> ``Lhx6``.
    Falls back to the raw title when the heuristic fails.
    """
    d = (describe or "").strip()
    if not d:
        return ""
    core = d.split("[")[0].strip()
    toks = re.findall(r"[A-Za-z0-9][A-Za-z0-9/_-]*", core)
    stop = {"protein", "factor", "isoform", "receptor", "kinase", "precursor",
            "variant", "like", "domain", "family", "member", "subunit"}
    # 1) symbol-like tokens (letters + digits, mixed case), e.g. Lhx6, TP63
    for t in toks:
        if re.fullmatch(r"[A-Za-z]{2,}[0-9]+[A-Za-z0-9]*", t):
            return t
    # 2) first alphabetic token that is not a generic descriptor
    for t in toks:
        low = t.lower()
        if t.isalpha() and len(t) >= 2 and low not in stop \
                and not low.startswith("isoform"):
            return t
    # 3) fall back to the first 2+ alnum token
    for t in toks:
        if len(t) >= 2 and not t.isdigit() and not t.startswith("("):
            return t
    return describe


def summarize_annotations(long_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse the long-format PlanMine parquet to one row per gene.

    Columns: gene_id_v6, contig_length, n_domains, dna_binding_domains,
    domains_all, go_terms, go_namespaces, blast_targets,
    planmine_human_ortholog_desc, rnai_markers_assayed (empty here).
    """
    rows = []
    bases = {}
    if "kind" in long_df.columns:
        for gid, sub in long_df[long_df["kind"] == "base"].groupby("gene_id_v6"):
            try:
                bases[gid] = int(sub.iloc[0]["contig_length"])
            except (ValueError, TypeError):
                bases[gid] = None

    for gid, sub in long_df.groupby("gene_id_v6"):
        go = sub[(sub["kind"] == "go")] if "kind" in sub else pd.DataFrame()
        dom = sub[(sub["kind"] == "domain")] if "kind" in sub else pd.DataFrame()
        blast = sub[(sub["kind"] == "blast")] if "kind" in sub else pd.DataFrame()
        go_names = [str(v) for v in go["value"].tolist() if str(v) not in ("nan", "")]
        go_ns = [str(v) for v in go["namespace"].tolist() if str(v) not in ("nan", "")]
        go_ids = [str(v) for v in go["key"].tolist() if str(v) not in ("nan", "")]
        dom_names = [str(v) for v in dom["value"].tolist() if str(v) not in ("nan", "")]
        dbs = sorted({d for d in dom_names if domain_short_name_is_dna_binding(d)})
        human_desc = ""
        for _, b in blast.iterrows():
            if "homo sapiens" in (b["value"] or "").lower():
                human_desc = str(b.get("description", ""))
                break
        rows.append({
            "gene_id_v6": gid,
            "contig_length": bases.get(gid),
            "n_domains": len(dom_names),
            "dna_binding_domains": ", ".join(dbs),
            "domains_all": "; ".join(sorted(set(dom_names))),
            "go_terms": "; ".join(go_names),
            "go_namespaces": "; ".join(go_ns),
            "go_ids": "; ".join(go_ids),
            "planmine_human_ortholog_desc": human_desc,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# scoring + selection
# ---------------------------------------------------------------------------
def prepare_candidates(rank: pd.DataFrame, mmc4: pd.DataFrame | None = None) -> pd.DataFrame:
    """Normalize a rank table and add King mmc4-derived features.

    Expected rank columns (from rank_neural.csv): gene_id, gene_name,
    integrated_score, proof_status, n_streams, rnai, ...
    Adds: better gene names, ``human_ortholog``, ``mmc4_tf_flag`` and empty
    slots for PlanMine annotation columns (filled later by merge_annotations).
    """
    r = rank.copy()
    r["gene_id"] = r["gene_id"].astype(str).str.strip()
    r["gene_name"] = r["gene_name"].fillna("").astype(str).str.strip()

    # King mmc4 catalog: better gene names + local human orthologs + TF flag
    r["human_ortholog"] = ""
    r["mmc4_tf_flag"] = ""
    r["gene_name_full"] = ""
    if mmc4 is not None and not mmc4.empty:
        m = mmc4.copy()
        m["Gene ID"] = m["Gene ID"].astype(str).str.strip()
        if "Planarian GenBank Gene Name" in m.columns:
            name_map = dict(zip(m["Gene ID"], m["Planarian GenBank Gene Name"].astype(str)))
            # keep a verbose copy; only fill a blank gene_name (never clobber the
            # clean pipeline symbol like "pitx" with a long FASTA-style title)
            r["gene_name_full"] = r["gene_id"].map(
                lambda g: name_map.get(g, "") if name_map.get(g, "") not in ("nan", "") else ""
            )
            r["gene_name"] = [
                (name_map.get(g, n) if (not str(n).strip() or str(n).strip() == "nan")
                 and name_map.get(g, "") not in ("nan", "") else n)
                for g, n in zip(r["gene_id"], r["gene_name"])
            ]
        if "Human Best Blast Hit" in m.columns:
            hu = dict(zip(m["Gene ID"], m["Human Best Blast Hit"].astype(str)))
            r["human_ortholog"] = r["gene_id"].map(
                lambda g: hu.get(g, "") if hu.get(g, "") != "nan" else ""
            )
        if "TF?" in m.columns:
            tfm = dict(zip(m["Gene ID"], m["TF?"].astype(str)))
            r["mmc4_tf_flag"] = r["gene_id"].map(
                lambda g: tfm.get(g, "") if tfm.get(g, "") != "nan" else ""
            )

    # PlanMine annotation columns (populated by a later merge; default empty)
    for col in (
        "dna_binding_domains", "domains_all", "go_terms", "go_namespaces",
        "planmine_human_ortholog_desc",
    ):
        if col not in r.columns:
            r[col] = ""
    if "rnai_phenotype_notes" not in r.columns:
        r["rnai_phenotype_notes"] = ""
    return r


def merge_annotations(score_df: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    """Join per-gene PlanMine annotation summaries onto candidate rows.

    Left-side default columns collide by name with the annotation columns
    (e.g. ``dna_binding_domains``); the annotation values win.
    """
    if ann is None or ann.empty:
        return score_df
    cols = [c for c in ann.columns if c != "gene_id_v6"]
    keep = score_df.copy()
    for c in cols:
        keep = keep.drop(columns=[c], errors="ignore")
    out = keep.merge(ann, left_on="gene_id", right_on="gene_id_v6", how="left")
    for c in cols:
        if c not in out.columns:
            out[c] = ""
        out[c] = out[c].fillna("")
    return out


def compute_composite(df: pd.DataFrame) -> pd.DataFrame:
    """Compute composite_score from the merged feature columns (mutating copy)."""
    out = df.copy()
    out["composite_score"] = out.apply(_composite_score, axis=1).round(4)
    return out


def _has_ortholog(row: pd.Series) -> bool:
    """True only when a real human-ortholog string is present.

    NaN (float) or the literal string "nan" must never count as evidence —
    pandas turns empty mmc4 cells into NaN, which is falsy for truthiness
    tests but truthy for `str(NaN)` -> "nan".
    """
    for col in ("human_ortholog", "planmine_human_ortholog_desc"):
        v = row.get(col, None)
        if pd.isna(v):
            continue
        s = str(v).strip()
        if s and s.lower() != "nan":
            return True
    return False


def _composite_score(row: pd.Series) -> float:
    base = 0.0
    try:
        base = float(row.get("integrated_score") or 0.0)
    except (TypeError, ValueError):
        base = 0.0
    bonus = 0.0
    doms = str(row.get("dna_binding_domains", "") or "").strip()
    if doms and doms.lower() not in ("nan", "none"):
        bonus += BONUS_TF_DOMAIN
    elif str(row.get("mmc4_tf_flag", "")).strip().upper() == "TF":
        bonus += BONUS_TF_DOMAIN
    has_neural_go = has_go_tf = False
    seen: set[str] = set()
    for term in str(row.get("go_terms", "") or "").split(";"):
        t = term.strip()
        if not t or t.lower() in ("nan", "none") or t.lower() in seen:
            continue
        seen.add(t.lower())
        is_neural, is_tf = go_term_flags(t)
        has_neural_go = has_neural_go or is_neural
        has_go_tf = has_go_tf or is_tf
    if has_neural_go:
        bonus += BONUS_GO_NEURAL
    if has_go_tf:
        bonus += BONUS_GO_TF
    if _has_ortholog(row):
        bonus += BONUS_HUMAN_ORTHOLOG
    if str(row.get("proof_status", "")).strip() == "known_rnai_validated":
        bonus += BONUS_BRAIN_PHENOTYPE
    return min(MAX_COMPOSITE, base + bonus)


def select_top(track_df: pd.DataFrame, n: int = 5) -> pd.DataFrame:
    """Take the top-n rows by composite_score (ties broken by the pipeline
    integrated_score, then n_streams).  Returns the rows ranked 1..n."""
    out = track_df.copy()
    for col in ("composite_score", "integrated_score"):
        if col not in out.columns:
            out[col] = 0.0
    out["composite_score"] = pd.to_numeric(
        out["composite_score"], errors="coerce"
    ).fillna(out["integrated_score"])
    out["integrated_score"] = pd.to_numeric(
        out["integrated_score"], errors="coerce"
    ).fillna(0.0)
    if "n_streams" not in out.columns:
        out["n_streams"] = 0
    out["n_streams"] = pd.to_numeric(
        out["n_streams"], errors="coerce"
    ).fillna(0).astype(int)
    out = out.sort_values(
        by=["composite_score", "integrated_score", "n_streams"],
        ascending=False,
    ).head(n).copy()
    out["rank"] = range(1, len(out) + 1)
    return out.reset_index(drop=True)


def assign_tracks(rank: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split into Track A (known_rnai_validated) and Track B (novel)."""
    a = rank[rank["proof_status"] == "known_rnai_validated"].copy()
    b = rank[rank["proof_status"] == "novel_candidate"].copy()
    return a, b


# ---------------------------------------------------------------------------
# RNAi screen context (King mmc5)
# ---------------------------------------------------------------------------
def rnai_marker_notes(
    mmc5: pd.DataFrame | None, gene_id: str, extras: str = ""
) -> str:
    """Build a human-readable RNAi phenotype note from the King mmc5 table.

    mmc5 rows list the silenced TF first and the cell-type markers assayed in
    the following columns.  Returns ``"screened; markers assayed: ..."``.
    """
    if mmc5 is None or mmc5.empty:
        return extras or "RNAi-validated (King 2024 mmc5); no marker detail available"
    tokens = set()
    for _, row in mmc5.iterrows():
        first = row.iloc[0]
        first_s = str(first).strip() if isinstance(first, str) else str(first)
        if first_s == gene_id or (gene_id in first_s):
            for cell in row.iloc[1:]:
                if isinstance(cell, str):
                    c = cell.strip()
                    if c and c != "nan":
                        tokens.add(c)
    prefix = extras if extras else "RNAi-validated (King 2024, mmc5)"
    if tokens:
        return f"{prefix}; phenotype screen markers assayed: {', '.join(sorted(tokens))}"
    return prefix