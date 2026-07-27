"""Gene identifier bridging between planarian genome builds.

The central engineering problem for NeuralTF is that Fincher et al. 2018
uses gene IDs from the **dd_Smed_v4** assembly while King 2024 and Plass
2018 use **dd_Smed_v6**. This module ships a :class:`BridgeTable` and the
helpers needed to load and build bridges from name-matching passes.

The framework never guesses identifier equivalence from numeric prefixes —
it requires an explicit bridge table.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd

from bioforge.core.logging import get_logger

logger = get_logger("evidence.gene_mapping")


@dataclass
class BridgeTable:
    """Mapping between gene identifiers across genome builds.

    The table is a thin wrapper around a :class:`pandas.DataFrame` with
    columns ``["gene_name", "v6_id", "v4_id"]``; any of ``v6_id``/``v4_id``
    may be NaN for genes only present in one build.
    """

    df: pd.DataFrame

    def __post_init__(self) -> None:
        required = {"gene_name", "v6_id", "v4_id"}
        missing = required - set(self.df.columns)
        if missing:
            raise ValueError(f"BridgeTable missing required columns: {sorted(missing)}")

    @property
    def n_rows(self) -> int:
        return int(len(self.df))

    @property
    def n_bridged(self) -> int:
        return int((self.df["v6_id"].notna() & self.df["v4_id"].notna()).sum())

    def v6_to_name(self, v6_id: str) -> Optional[str]:
        sub = self.df.loc[self.df["v6_id"] == v6_id, "gene_name"]
        if sub.empty:
            return None
        return str(sub.iloc[0])

    def v4_to_name(self, v4_id: str) -> Optional[str]:
        sub = self.df.loc[self.df["v4_id"] == v4_id, "gene_name"]
        if sub.empty:
            return None
        return str(sub.iloc[0])

    def v4_to_v6(self, v4_id: str) -> Optional[str]:
        sub = self.df.loc[self.df["v4_id"] == v4_id, "v6_id"]
        if sub.empty or pd.isna(sub.iloc[0]):
            return None
        return str(sub.iloc[0])

    def v6_to_v4(self, v6_id: str) -> Optional[str]:
        sub = self.df.loc[self.df["v6_id"] == v6_id, "v4_id"]
        if sub.empty or pd.isna(sub.iloc[0]):
            return None
        return str(sub.iloc[0])


def load_bridge(path: str | Path) -> BridgeTable:
    """Load a bridge table from CSV (``gene_name,v6_id,v4_id``)."""
    df = pd.read_csv(path)
    return BridgeTable(df=df)


def build_bridge_from_names(
    v6_table: pd.DataFrame,
    v4_table: pd.DataFrame,
    *,
    v6_name_col: str = "gene_name",
    v6_id_col: str = "v6_id",
    v4_name_col: str = "gene_name",
    v4_id_col: str = "v4_id",
) -> BridgeTable:
    """Outer-join two name-id tables on ``gene_name``.

    Parameters
    ----------
    v6_table
        DataFrame mapping ``gene_name`` → ``v6_id`` for the v6 build.
    v4_table
        DataFrame mapping ``gene_name`` → ``v4_id`` for the v4 build.

    Notes
    -----
    Names must be unique within each input table; duplicate names are
    dropped with a warning so the bridge is deterministic.
    """
    for name, df, name_col in (
        ("v6", v6_table, v6_name_col),
        ("v4", v4_table, v4_name_col),
    ):
        if df[name_col].duplicated().any():
            dups = df.loc[df[name_col].duplicated(), name_col].tolist()
            logger.warning("dropping duplicate %s gene names: %s", name, dups[:5])
    v6 = v6_table[[v6_name_col, v6_id_col]].rename(
        columns={v6_name_col: "gene_name", v6_id_col: "v6_id"}
    )
    v6 = v6.drop_duplicates(subset="gene_name", keep=False)
    v4 = v4_table[[v4_name_col, v4_id_col]].rename(
        columns={v4_name_col: "gene_name", v4_id_col: "v4_id"}
    )
    v4 = v4.drop_duplicates(subset="gene_name", keep=False)
    merged = v6.merge(v4, on="gene_name", how="outer")
    logger.info(
        "built bridge: %d rows, %d fully bridged",
        len(merged),
        (merged["v6_id"].notna() & merged["v4_id"].notna()).sum(),
    )
    return BridgeTable(df=merged)
