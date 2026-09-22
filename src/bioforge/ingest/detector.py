"""Content-sniffing format detection.

Decision is based on a short peek at file bytes plus path metadata. We
try sniffers in priority order; the first that returns a confident answer
wins. Unreadable inputs produce :class:`UnknownFormatError` rather than
a partial / wrong result.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Union

from bioforge.core.logging import get_logger
from bioforge.ingest.readers import (
    read_10x_mtx,
    read_csv_matrix,
    read_dge_gz,
    read_h5ad,
    read_tsv_matrix,
)

logger = get_logger("ingest.detector")


HDF5_MAGIC = b"\x89HDF\r\n\x1a\n"
GZIP_MAGIC = b"\x1f\x8b"


@dataclass
class FormatInfo:
    format: str  # "h5ad" | "10x_mtx" | "dge_gz" | "tsv" | "csv" | "csv_gz"
    reader: Callable[..., object]
    needs_directory: bool = False


class UnknownFormatError(ValueError):
    """Raised when no reader accepts the supplied dataset."""


def _gz_csv_reader(path):
    import pandas as pd
    from bioforge.ingest.readers import _df_genes_to_anndata
    df = pd.read_csv(path, index_col=0, compression="gzip")
    return _df_genes_to_anndata(df, path)


class FormatDetector:
    """Detect dataset format by sniffing bytes and inspecting the path."""

    def detect(self, path: Union[str, Path]) -> FormatInfo:
        p = Path(path)
        if p.is_dir():
            if (p / "matrix.mtx").exists() or (p / "matrix.mtx.gz").exists():
                return FormatInfo("10x_mtx", read_10x_mtx, needs_directory=True)
            raise UnknownFormatError(
                f"directory {p} doesn't look like a 10x mtx layout"
            )
        if not p.exists():
            raise UnknownFormatError(f"file not found: {p}")

        suffixes = "".join(p.suffixes[-2:]).lower()
        if suffixes.endswith(".h5ad") or suffixes.endswith(".h5"):
            return FormatInfo("h5ad", read_h5ad)
        if suffixes.endswith(".mtx.gz") or suffixes.endswith(".mtx"):
            raise UnknownFormatError(
                f"matrix.mtx should be inside a 10x-style directory; got {p.name}"
            )

        head = _peek(p, n=512)

        if suffixes.endswith(".gz"):
            # _peek transparently gunzips a small chunk for sniffing, so
            # if the underlying compressed content is actually a gzipped HDF5
            # we wouldn't be here (`.h5ad.gz` isn't a BioForge standard).
            # Treat any text-looking content as a zipped matrix.
            return self._detect_text_gz(head, suffixes)

        if suffixes.endswith(".csv"):
            return FormatInfo("csv", read_csv_matrix)
        if suffixes.endswith(".tsv") or suffixes.endswith(".txt"):
            return FormatInfo("tsv", read_tsv_matrix)

        # Unknown extension — sniff by content.
        if head.startswith(HDF5_MAGIC):
            return FormatInfo("h5ad", read_h5ad)
        if head.startswith(GZIP_MAGIC):
            return self._detect_text_gz(head, suffixes)

        return self._detect_text_plain(head)

    def _detect_text_gz(self, head: bytes, suffixes: str) -> FormatInfo:
        text = head.decode("utf-8", "replace")
        first_line = text.splitlines()[0] if text else ""
        n_tabs = first_line.count("\t")
        n_commas = first_line.count(",")
        if n_tabs > 0 and n_tabs >= n_commas:
            return FormatInfo("dge_gz", read_dge_gz)
        if n_commas > 0:
            logger.warning("gzipped CSV detected; using csv+gzip adapter")
            return FormatInfo("csv_gz", _gz_csv_reader)
        raise UnknownFormatError("gzip text but no separator found in header")

    def _detect_text_plain(self, head: bytes) -> FormatInfo:
        text = head.decode("utf-8", "replace")
        first_line = text.splitlines()[0] if text else ""
        n_tabs = first_line.count("\t")
        n_commas = first_line.count(",")
        if n_tabs > 0 and n_tabs >= n_commas:
            return FormatInfo("tsv", read_tsv_matrix)
        if n_commas > 0:
            return FormatInfo("csv", read_csv_matrix)
        raise UnknownFormatError(
            "file looks like text but has neither tabs nor commas in its header"
        )


def _peek(path: Path, n: int = 512) -> bytes:
    with open(path, "rb") as fh:
        head = fh.read(n)
    if head.startswith(GZIP_MAGIC):
        with gzip.open(path, "rb") as gf:
            head = gf.read(n)
    return head
