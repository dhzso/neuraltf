"""BioForge dataset ingestion (Layer 8C).

Single entry point :func:`ingest_dataset` accepts GEO/SRA accessions, URLs,
or local paths and returns an :class:`anndata.AnnData` -- auto-detecting
format from content sniffing (not just extension). Unreadable formats
degrade gracefully and surface as :class:`UnknownFormatError` so the
caller (UI/CLI) can render a helpful message instead of a stack trace.

Optional FASTQ-to-matrix orchestration lives in :mod:`bioforge.ingest.fastq`
and is loaded lazily so the heavy kb-python/salmon toolchain doesn't need
to be installed unless the user opts in.
"""
from bioforge.ingest.api import (
    DatasetSource,
    UnknownFormatError,
    ingest_dataset,
    resolve_source,
)
from bioforge.ingest.detector import FormatDetector, FormatInfo
from bioforge.ingest.readers import (
    read_csv_matrix,
    read_dge_gz,
    read_h5ad,
    read_10x_mtx,
    read_tsv_matrix,
)
from bioforge.ingest.resolver import (
    is_accession,
    is_url,
)

__all__ = [
    "DatasetSource",
    "UnknownFormatError",
    "ingest_dataset",
    "FormatInfo",
    "FormatDetector",
    "read_dge_gz",
    "read_h5ad",
    "read_10x_mtx",
    "read_csv_matrix",
    "read_tsv_matrix",
    "is_accession",
    "is_url",
    "resolve_source",
]
