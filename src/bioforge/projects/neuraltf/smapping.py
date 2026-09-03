"""Unified gene ID mapping across SMED, h1SMcG, and dd_Smed_v6 systems.

Provides lazy-loaded mappings between the three gene identifier conventions
used by Cui 2023 (SMED), Perez 2025 (h1SMcG), and the NeuralTF pipeline
(dd_Smed_v6).
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd

# Resolve relative to this file's location, not CWD
_REPO = Path(__file__).resolve().parents[4]
_RAW = _REPO / "datasets" / "raw"
_ROSETTA = _RAW / "smed_20140614.mapping.rosettastone.2020.txt"
_MOESM5 = (
    _RAW
    / "Supplementary_Data_ Perez_2025"
    / "41467_2025_65712_MOESM5_ESM.xlsx"
)


@lru_cache(maxsize=1)
def _load_rosetta() -> pd.DataFrame:
    """Load Rosetta Stone and keep only SMED -> dd_Smed_v6 rows."""
    df = pd.read_csv(
        _ROSETTA,
        sep="\t",
        header=None,
        names=["ref_id", "seq_id", "transcriptome_id"],
        usecols=[0, 1, 2],
        dtype=str,
        on_bad_lines="skip",
    )
    v6 = df[df["transcriptome_id"].str.contains("v6", na=False)].copy()
    v6["smed_id"] = v6["ref_id"].str.strip()
    v6["v6_id"] = v6["seq_id"].str.strip()
    return v6[["smed_id", "v6_id"]].dropna().drop_duplicates()


@lru_cache(maxsize=1)
def _load_moism5() -> pd.DataFrame:
    """Load MOESM5 and extract h1SMcG <-> dd_Smed_v6 mappings."""
    raw = pd.read_excel(_MOESM5, sheet_name=0, dtype=str)
    cols = raw.columns.tolist()
    gene_col = cols[0]
    rbh_col = next((c for c in cols if "1:1" in c and "v6" in c.lower()), None)
    collapsed_col = next((c for c in cols if "Similar" in c and "v6" in c.lower()), None)
    tf_class_col = next((c for c in cols if "TF Class" in c and "Perez" in c), None)
    tf_module_col = next((c for c in cols if "module of TF" in c and "Perez" in c), None)

    records = raw.to_dict(orient="records")
    rows = []
    for r in records:
        h1 = str(r.get(gene_col, "")).strip()
        if not h1 or h1 == "nan":
            continue
        v6_rbh = str(r.get(rbh_col, "")).strip() if rbh_col else ""
        v6_coll = str(r.get(collapsed_col, "")).strip() if collapsed_col else ""
        tf_cls = str(r.get(tf_class_col, "")).strip() if tf_class_col else ""
        tf_mod = str(r.get(tf_module_col, "")).strip() if tf_module_col else ""
        v6_ids = set()
        for v in [v6_rbh, v6_coll]:
            if v and v != "nan":
                for part in v.replace(";", ",").split(","):
                    p = part.strip()
                    if p.startswith("dd_Smed_v6"):
                        v6_ids.add(p)
        rows.append({
            "h1smcg_id": h1,
            "v6_rbh": v6_rbh if v6_rbh != "nan" else "",
            "v6_all": sorted(v6_ids),
            "tf_class": tf_cls if tf_cls != "nan" else "",
            "tf_module": tf_mod if tf_mod != "nan" else "",
        })
    return pd.DataFrame(rows)


@lru_cache(maxsize=1)
def _get_v6_maps() -> tuple[dict[str, str], dict[str, str]]:
    """Build fast O(1) hash maps for v6 -> h1SMcG and h1SMcG -> TF Class."""
    df = _load_moism5()
    v6_to_h1: dict[str, str] = {}
    h1_to_class: dict[str, str] = {}
    for r in df.itertuples(index=False):
        h1 = r.h1smcg_id
        if r.tf_class:
            h1_to_class[h1] = r.tf_class
        if r.v6_rbh:
            v6_to_h1[r.v6_rbh] = h1
        for v in r.v6_all:
            if v not in v6_to_h1:
                v6_to_h1[v] = h1
    return v6_to_h1, h1_to_class


@lru_cache(maxsize=1)
def _get_v6_maps_rbh() -> dict[str, str]:
    """RBH-restricted v6 -> h1SMcG map.

    ~14.4k of 25k v6 IDs are claimed by more than one h1SMcG via the
    collapsed (`Similar`) column; first-wins picks are arbitrary there.
    The 1:1 reciprocal-best-hit column is unambiguous, so influence and
    lineage attribution should use only these rows.
    """
    df = _load_moism5()
    v6_to_h1: dict[str, str] = {}
    for r in df.itertuples(index=False):
        if r.v6_rbh and r.v6_rbh not in v6_to_h1:
            v6_to_h1[r.v6_rbh] = r.h1smcg_id
    return v6_to_h1


def smed_to_v6(smed_id: str) -> list[str]:
    """Map a SMED300* ID to a list of dd_Smed_v6 IDs."""
    df = _load_rosetta()
    hits = df[df["smed_id"] == smed_id.strip()]["v6_id"].tolist()
    return sorted(set(hits))


def v6_to_smed(v6_id: str) -> list[str]:
    """Map a dd_Smed_v6 ID to a list of SMED300* IDs."""
    df = _load_rosetta()
    hits = df[df["v6_id"] == v6_id.strip()]["smed_id"].tolist()
    return sorted(set(hits))


def h1smcg_to_v6(h1smcg_id: str) -> list[str]:
    """Map a h1SMcG ID to a list of dd_Smed_v6 IDs."""
    df = _load_moism5()
    row = df[df["h1smcg_id"] == h1smcg_id.strip()]
    if row.empty:
        return []
    return row.iloc[0]["v6_all"]


def v6_to_h1smcg(v6_id: str) -> str | None:
    """Map a dd_Smed_v6 ID to its h1SMcG ID in O(1) time."""
    v6_to_h1, _ = _get_v6_maps()
    return v6_to_h1.get(v6_id.strip())


def v6_to_perez_tf_class(v6_id: str) -> str | None:
    """Look up the Perez 2025 TF class for a dd_Smed_v6 gene in O(1) time."""
    v6_to_h1, h1_to_class = _get_v6_maps()
    h1 = v6_to_h1.get(v6_id.strip())
    return h1_to_class.get(h1) if h1 else None


def is_perez_tf(v6_id: str) -> bool:
    """Check if a gene is classified as a TF by Perez 2025."""
    cls = v6_to_perez_tf_class(v6_id)
    return bool(cls and cls.lower() not in ("", "nan", "none"))


def batch_v6_to_h1smcg(v6_ids: list[str], *, rbh_only: bool = False) -> dict[str, str | None]:
    """Map a list of dd_Smed_v6 IDs to h1SMcG IDs in O(1) time.

    Set ``rbh_only=True`` to restrict to unambiguous 1:1 reciprocal best
    hits (recommended for scoring evidence streams).
    """
    if rbh_only:
        v6_to_h1 = _get_v6_maps_rbh()
    else:
        v6_to_h1, _ = _get_v6_maps()
    return {v6: v6_to_h1.get(v6) for v6 in v6_ids}


def mapping_stats() -> dict:
    """Return mapping quality statistics for all three maps.

    Returns
    -------
    dict with keys:
        rosetta: {total_smed, total_v6, mapped_smed, mapped_v6, rate_smed_to_v6, rate_v6_to_smed, one_to_many}
        moesm5: {total_h1smcg, total_v6_rbh, total_v6_all, mapped_h1smcg, mapped_v6, rate_h1smcg_to_v6, rate_v6_to_h1smcg}
    """
    rosetta = _load_rosetta()
    moesm5 = _load_moism5()

    # Rosetta Stone stats
    total_smed = rosetta["smed_id"].nunique()
    total_v6 = rosetta["v6_id"].nunique()
    mapped_smed = total_smed  # all SMED in table map to at least one v6
    mapped_v6 = total_v6      # all v6 in table map to at least one SMED

    # One-to-many counts
    smed_to_v6_counts = rosetta.groupby("smed_id")["v6_id"].nunique()
    v6_to_smed_counts = rosetta.groupby("v6_id")["smed_id"].nunique()
    smed_one_to_many = (smed_to_v6_counts > 1).sum()
    v6_one_to_many = (v6_to_smed_counts > 1).sum()

    # MOESM5 stats
    total_h1smcg = moesm5["h1smcg_id"].nunique()
    total_v6_rbh = moesm5["v6_rbh"].replace("", pd.NA).dropna().nunique()
    all_v6 = set()
    for v6_list in moesm5["v6_all"]:
        all_v6.update(v6_list)
    total_v6_all = len(all_v6)

    mapped_h1smcg = moesm5[moesm5["v6_all"].apply(len) > 0]["h1smcg_id"].nunique()
    mapped_v6_moesm5 = len([v for v in all_v6 if v6_to_h1smcg(v) is not None])

    # One-to-many in MOESM5
    h1_to_v6_counts = moesm5["v6_all"].apply(len)
    h1_one_to_many = (h1_to_v6_counts > 1).sum()

    return {
        "rosetta": {
            "total_smed": int(total_smed),
            "total_v6": int(total_v6),
            "mapped_smed": int(mapped_smed),
            "mapped_v6": int(mapped_v6),
            "rate_smed_to_v6": 1.0,
            "rate_v6_to_smed": 1.0,
            "smed_one_to_many": int(smed_one_to_many),
            "v6_one_to_many": int(v6_one_to_many),
        },
        "moesm5": {
            "total_h1smcg": int(total_h1smcg),
            "total_v6_rbh": int(total_v6_rbh),
            "total_v6_all": int(total_v6_all),
            "mapped_h1smcg": int(mapped_h1smcg),
            "mapped_v6": int(mapped_v6_moesm5),
            "rate_h1smcg_to_v6": float(mapped_h1smcg / total_h1smcg) if total_h1smcg > 0 else 0.0,
            "rate_v6_to_h1smcg": float(mapped_v6_moesm5 / total_v6_all) if total_v6_all > 0 else 0.0,
            "h1_one_to_many": int(h1_one_to_many),
        },
    }
