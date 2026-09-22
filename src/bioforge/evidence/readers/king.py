"""Reader for King 2024 supplementary xlsx files.

The King supplement ships a single `mmc4.xlsx` TF catalog, an `mmc5.xlsx`
RNAi phenotype table, an `mmc6.xlsx` neural TF-pair correlation table, and
the `mmc7.xlsx` TF atlas (enrichment + p-values + log2FC for both
post-mitotic G0 progenitors and X1 neoblasts).

The xlsx files have multi-line descriptive headers above the real column
row, so readers must skip a fixed number of leading rows per sheet. We
hard-code those skips based on inspection of the supplied files
(``datasets/raw/Supplementary_Data_ King_2024``).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from bioforge.core.logging import get_logger

logger = get_logger("evidence.readers.king")


def _read_sheet(path: str | Path, sheet: str, header_row: int) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=sheet, header=header_row)
    df = df.dropna(how="all").reset_index(drop=True)
    return df


# ---- mmc4: TF catalog ----------------------------------------------------


def read_tf_catalog(path: str | Path) -> pd.DataFrame:
    """Read the King mmc4 TF catalog.

    Returns
    -------
    DataFrame with columns:
        ``gene_id`` (dd_Smed_v6), ``human_hit_acc``, ``human_hit``,
        ``e_value``, ``planarian_name``, ``percent_identity``,
        ``is_tf`` (bool), ``is_prior_fstf`` (bool)
    """
    df = _read_sheet(path, sheet="TF", header_row=0)
    out = df.rename(columns={
        "Gene ID": "gene_id",
        "Human Best Blast Hit Accession Number": "human_hit_acc",
        "Human Best Blast Hit": "human_hit",
        "E-value": "e_value",
        "Planarian GenBank Gene Name": "planarian_name",
        "Percent Nucloetide Identity": "percent_identity",
        "TF?": "is_tf",
        "FSTF?": "is_prior_fstf",
    })
    out["is_tf"] = out["is_tf"].astype(str).str.strip().str.upper().eq("TF")
    out["is_prior_fstf"] = (
        out["is_prior_fstf"].astype(str).str.strip().str.lower().eq("yes")
    )
    logger.info("loaded TF catalog: %d entries", len(out))
    return out


# ---- mmc5: RNAi phenotypes -----------------------------------------------


def read_rnai_table(path: str | Path) -> pd.DataFrame:
    """Read the King mmc5 RNAi phenotype table.

    The first data column (col 0) holds the knocked-down gene ID; remaining
    columns (often named after marker genes or "Cell-type markers tested")
    list markers that were tested for a phenotype. The reader returns a
    long-format DataFrame with columns ``fstf_rnai`` and ``marker``.
    """
    raw = pd.read_excel(path, sheet_name=0, header=None)
    # Rows 0-3 are headings; data starts at row 4 (0-indexed).
    data = raw.iloc[4:].reset_index(drop=True)
    records = []
    for _, row in data.iterrows():
        gene = row.iloc[0]
        if pd.isna(gene):
            continue
        # Remaining cells in the row contain marker IDs; pd.isna guards
        # against empty/None cells (some rows are ragged).
        for cell in row.iloc[1:]:
            if pd.notna(cell) and isinstance(cell, str) and cell.strip():
                records.append({"fstf_rnai": str(gene), "marker": cell.strip()})
    out = pd.DataFrame(records, columns=["fstf_rnai", "marker"])
    logger.info("loaded RNAi table: %d (FSTF, marker) pairs", len(out))
    return out


# ---- mmc6: neural TF pair correlations -----------------------------------


def read_correlations(path: str | Path) -> pd.DataFrame:
    """Read the King mmc6 neural TF pair correlations.

    The header row is row 3 (0-indexed). Columns: ``TF1``, ``TF2``,
    ``X1 Correlation``, ``G0 Correlation``, ``G0 Cluster``.
    """
    df = _read_sheet(path, sheet="Sheet1", header_row=3)
    df = df.rename(columns={
        "TF1": "tf1",
        "TF2": "tf2",
        "X1 Correlation": "x1_corr",
        "G0 Correlation": "g0_corr",
        "G0 Cluster": "g0_cluster",
    })
    logger.info("loaded TF pair correlations: %d pairs", len(df))
    return df


# ---- mmc7: TF atlas ------------------------------------------------------


_MMC7_SHEETS = {
    "g0_atlas": ("G0 Progenitor TF Atlas", 4),
    "g0_pvalues": ("G0 Progenitor Pvalues", 3),
    "g0_log2fc": ("G0 Progenitor Log2FC", 3),
    "x1_atlas": ("X1 TF Atlas", 4),
    "x1_pvalues": ("X1 Pvalues", 3),
    "x1_log2fc": ("X1 Log2FC", 3),
}


def read_king_atlas(path: str | Path) -> dict[str, pd.DataFrame]:
    """Read the six core sheets of mmc7 (G0 and X1 atlas+values).

    Returns a dict keyed by the short names in ``_MMC7_SHEETS``; each value
    is the raw sheet content with its first data column renamed to
    ``subcluster`` when present.
    """
    out: dict[str, pd.DataFrame] = {}
    for short, (sheet, header_row) in _MMC7_SHEETS.items():
        df = _read_sheet(path, sheet=sheet, header_row=header_row)
        # Standardize the first column name regardless of its xlsx label.
        first_col = df.columns[0]
        df = df.rename(columns={first_col: "subcluster"})
        out[short] = df
        logger.info("mmc7 sheet '%s' -> %d rows", short, len(df))
    return out
