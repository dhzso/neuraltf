"""Top-level dataset ingestion API.

This is the public entry point: pass a GEO/SRA accession, a URL, or a local
path; get back an :class:`anndata.AnnData`. Unknown formats raise
:class:`UnknownFormatError`, which the UI catches and renders as a user-
friendly message rather than a stack trace.
"""
from __future__ import annotations

import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

import anndata as ad

from bioforge.core.logging import get_logger
from bioforge.ingest.detector import FormatDetector, FormatInfo, UnknownFormatError
from bioforge.ingest.resolver import classify_source, is_accession, is_url

logger = get_logger("ingest.api")


@dataclass
class DatasetSource:
    raw: str
    kind: str  # "accession" | "url" | "local_path"
    local_path: Optional[Path] = None
    detection_note: str = ""


def resolve_source(source: str, *, dest_dir: Optional[Union[str, Path]] = None) -> DatasetSource:
    """Return a :class:`DatasetSource` describing what's pointed at by `source`.

    For accessions and URLs we don't fetch yet (no network); for local
    paths we record the resolved Path. A workflow step or the UI can
    then call :func:`ingest_dataset` to do format detection and reading.
    """
    kind = classify_source(source)
    if kind == "local_path":
        return DatasetSource(raw=source, kind=kind, local_path=Path(source).expanduser())
    return DatasetSource(raw=source, kind=kind, detection_note="deferred fetch")


_DOWNLOADABLE_EXTENSIONS = (".h5ad", ".csv", ".tsv", ".txt", ".txt.gz",
                            ".csv.gz", ".tsv.gz", ".h5")


def _ensure_local(src: DatasetSource, dest_dir: Path) -> Path:
    """For url sources, download to dest_dir; for accession, raise (handled
    by an explicit later step that hits GEO/SRA APIs)."""
    if src.local_path is not None:
        return src.local_path
    if src.kind == "url":
        dest_dir.mkdir(parents=True, exist_ok=True)
        name = src.raw.rsplit("/", 1)[-1] or "dataset.bin"
        # If the URL has a query string, strip it.
        name = name.split("?", 1)[0]
        dest = dest_dir / name
        logger.info("downloading %s -> %s", src.raw, dest)
        req = urllib.request.Request(src.raw, headers={"User-Agent": "BioForge/0.1"})
        with urllib.request.urlopen(req) as resp, open(dest, "wb") as fh:
            fh.write(resp.read())
        return dest
    if src.kind == "accession":
        # First implementation: pretend the accession maps to nothing locally;
        # throwing here lets the UI show a friendly "download not implemented
        # for this accession; please use `geo`/`sra-tools` to fetch first"
        # message. This is the graceful degraded path promised in ADR-0003.
        raise UnknownFormatError(
            f"Direct SRA/GEO download for accession '{src.raw}' is not yet "
            "implemented. Please run `geo` (NCBI Gene Expression Omnibus "
            "`wget` helper) or `sra-toolkit fasterq-dump` first and pass the "
            "local file path to BioForge."
        )
    raise UnknownFormatError(f"cannot resolve source of kind {src.kind!r}")


def ingest_dataset(
    source: Union[str, DatasetSource],
    *,
    dest_dir: Optional[Union[str, Path]] = None,
    detector: Optional[FormatDetector] = None,
) -> ad.AnnData:
    """Ingest a dataset and return an :class:`anndata.AnnData`.

    The function gracefully fails on unknown formats; the caller should
    catch :class:`UnknownFormatError` and surface a helpful UI/CLI message.
    """
    if isinstance(source, str):
        src = resolve_source(source)
    else:
        src = source
    dest = Path(dest_dir) if dest_dir else Path.cwd() / "datasets" / "cache"
    path = _ensure_local(src, dest)
    detector = detector or FormatDetector()
    info: FormatInfo = detector.detect(path)
    logger.info("detected format '%s' for %s", info.format, path)
    adata = info.reader(path)
    adata.uns.setdefault("bioforge_source", src.raw)
    adata.uns.setdefault("bioforge_format", info.format)
    return adata
