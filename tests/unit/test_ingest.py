"""Unit tests for bioforge.ingest (Layer 8C)."""
from __future__ import annotations

import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from bioforge.ingest import (
    FormatDetector,
    UnknownFormatError,
    ingest_dataset,
    is_accession,
    is_url,
    read_csv_matrix,
    read_dge_gz,
    read_h5ad,
    read_10x_mtx,
    read_tsv_matrix,
    resolve_source,
)


# ---------------------------------------------------------------------------
# Source type heuristics
# ---------------------------------------------------------------------------


def test_is_accession_matches_geo_and_sra() -> None:
    for s in ("GSE12345", "GSM22", "SRP3", "SRR9", "SRX11"):
        assert is_accession(s)
    assert not is_accession("not_an_accession")
    assert not is_accession("")


def test_is_url_recognizes_http_https_ftp() -> None:
    assert is_url("https://example.org/x.h5ad")
    assert is_url("ftp://ftp.ncbi.nlm.nih.gov/x.txt.gz")
    assert not is_url("GSE12345")


def test_resolve_source_for_local_path() -> None:
    src = resolve_source("datasets/raw/x.h5ad")
    assert src.kind == "local_path"
    assert src.local_path is not None


def test_resolve_source_for_accession_defers_fetch() -> None:
    src = resolve_source("GSE12345")
    assert src.kind == "accession"


# ---------------------------------------------------------------------------
# FormatDetector — content sniffing
# ---------------------------------------------------------------------------


def test_detector_for_h5ad_extension(tmp_path: Path) -> None:
    p = tmp_path / "x.h5ad"
    p.write_bytes(b"\x89HDF\r\n\x1a\n" + b"\x00" * 16)
    info = FormatDetector().detect(p)
    assert info.format == "h5ad"


def test_detector_for_dge_gz(tmp_path: Path) -> None:
    df = pd.DataFrame({"cell_1": [1, 2], "cell_2": [3, 4]},
                      index=["g1", "g2"])
    p = tmp_path / "dge.txt.gz"
    with gzip.open(p, "wb") as fh:
        fh.write(b"gene\tcell_1\tcell_2\n")
        for gene, row in df.iterrows():
            fh.write(f"{gene}\t{row['cell_1']}\t{row['cell_2']}\n".encode())
    info = FormatDetector().detect(p)
    assert info.format == "dge_gz"


def test_detector_for_csv_without_extension(tmp_path: Path) -> None:
    p = tmp_path / "noext"
    p.write_text("gene,c1,c2\ng1,1,2\ng2,3,4\n")
    info = FormatDetector().detect(p)
    assert info.format == "csv"


def test_detector_for_tsv_without_extension(tmp_path: Path) -> None:
    p = tmp_path / "noext"
    p.write_text("gene\tc1\tc2\ng1\t1\t2\ng2\t3\t4\n")
    info = FormatDetector().detect(p)
    assert info.format == "tsv"


def test_detector_for_10x_mtx_directory(tmp_path: Path) -> None:
    d = tmp_path / "10x"
    d.mkdir()
    (d / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate real general\n")
    (d / "barcodes.tsv").write_text("AAAC\n")
    (d / "features.tsv").write_text("g1\tg1\tgene\n")
    info = FormatDetector().detect(d)
    assert info.format == "10x_mtx"
    assert info.needs_directory is True


def test_detector_raises_on_missing_file(tmp_path: Path) -> None:
    with pytest.raises(UnknownFormatError):
        FormatDetector().detect(tmp_path / "nonexistent.h5ad")


def test_detector_raises_on_unsupported_extension_text(tmp_path: Path) -> None:
    p = tmp_path / "weird.unknown цифров.bin"
    p.write_bytes(b"\x00\x01\x02")
    with pytest.raises(UnknownFormatError):
        FormatDetector().detect(p)


def test_detector_raises_on_loose_mtx_file(tmp_path: Path) -> None:
    p = tmp_path / "matrix.mtx"
    p.write_text("whatever")
    with pytest.raises(UnknownFormatError) as exc:
        FormatDetector().detect(p)
    assert "10x-style directory" in str(exc.value)


# ---------------------------------------------------------------------------
# Readers — sanity round-trips
# ---------------------------------------------------------------------------


def test_read_csv_matrix_round_trip(tmp_path: Path) -> None:
    df = pd.DataFrame({"c1": [1, 2], "c2": [3, 4]}, index=["g1", "g2"])
    p = tmp_path / "m.csv"
    df.to_csv(p)
    adata = read_csv_matrix(p)
    assert adata.n_obs == 2 and adata.n_vars == 2
    assert list(adata.var_names) == ["g1", "g2"]


def test_read_dge_gz_round_trip(tmp_path: Path) -> None:
    df = pd.DataFrame({"c1": [1, 2], "c2": [3, 4]}, index=["g1", "g2"])
    p = tmp_path / "dge.txt.gz"
    with gzip.open(p, "wb") as fh:
        df.to_csv(fh, sep="\t")
    adata = read_dge_gz(p)
    assert adata.n_obs == 2 and adata.n_vars == 2


def test_read_h5ad_round_trip(tmp_path: Path) -> None:
    import anndata as _ad
    import pandas as _pd
    _ad.settings.allow_write_nullable_strings = True
    _pd.set_option("mode.string_storage", "python")
    obs = _pd.DataFrame(index=np.array([f"c{i}" for i in range(3)], dtype=object))
    var = _pd.DataFrame(index=np.array([f"g{i}" for i in range(3)], dtype=object))
    adata = ad.AnnData(
        X=np.eye(3, dtype=np.float32),
        obs=obs,
        var=var,
    )
    p = tmp_path / "x.h5ad"
    adata.write_h5ad(p)
    out = read_h5ad(p)
    assert out.n_obs == 3


# ---------------------------------------------------------------------------
# ingest_dataset end to end (local file only — no network)
# ---------------------------------------------------------------------------


def test_ingest_dataset_local_csv(tmp_path: Path) -> None:
    df = pd.DataFrame({"c1": [1, 2], "c2": [3, 4]}, index=["g1", "g2"])
    p = tmp_path / "m.csv"
    df.to_csv(p)
    out = ingest_dataset(str(p), dest_dir=tmp_path)
    assert isinstance(out, ad.AnnData)
    assert out.n_obs == 2
    assert out.uns["bioforge_format"] == "csv"


def test_ingest_dataset_unknown_format_raises_friendly(tmp_path: Path) -> None:
    p = tmp_path / "weird.bin"
    p.write_bytes(b"\x00\x01\x02")
    with pytest.raises(UnknownFormatError):
        ingest_dataset(str(p), dest_dir=tmp_path)


def test_ingest_dataset_accession_graceful_error(tmp_path: Path) -> None:
    # First iteration: direct accession fetch not implemented; surface a
    # user-friendly error rather than crashing.
    with pytest.raises(UnknownFormatError) as exc:
        ingest_dataset("GSE12345", dest_dir=tmp_path)
    assert "geo" in str(exc.value).lower() or "sra" in str(exc.value).lower()
