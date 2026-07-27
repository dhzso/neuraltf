"""Unit tests for bioforge.evidence.readers (King/Fincher/Plass loaders).

We build synthetic xlsx / txt files in tmp_path rather than depending on
the (git-ignored) raw datasets so the test suite is hermetic.
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import openpyxl
import pandas as pd
import pytest

from bioforge.evidence.readers import fincher, king, plass


# ---------------------------------------------------------------------------
# helpers — build a minimal King-style xlsx in memory
# ---------------------------------------------------------------------------


def _make_empty_workbook() -> openpyxl.Workbook:
    return openpyxl.Workbook()


def _write_tf_catalog(sh: openpyxl.worksheet.worksheet.Worksheet) -> None:
    # Header row (row 1, 1-indexed) matches King mmc4 sheet "TF"
    headers = [
        "Gene ID", "Human Best Blast Hit Accession Number", "Human Best Blast Hit",
        "E-value", "Planarian GenBank Gene Name", "Percent Nucloetide Identity",
        "TF?", "FSTF?",
    ]
    sh.append(headers)
    sh.append([
        "dd_Smed_v6_21801_0_1", "ref|XP_011537084.1|",
        "ALX homeobox protein 1 isoform X1 [Homo sapiens]", 0.0,
        "Schmidtea mediterranea aristaless-like homeobox transcription factor (arx) mRNA, complete cds",
        99.932, "TF", "Yes",
    ])
    sh.append([
        "dd_Smed_v6_1399_0_1", None, None, None,
        "hypothetical protein mRNA", 99.444, "NoBB", "No",
    ])


def _write_rnai(sh: openpyxl.worksheet.worksheet.Worksheet) -> None:
    # Rows 1-3 are descriptive header text, row 4 is the column header,
    # data starts row 5 (1-indexed).
    for _ in range(3):
        sh.append([" Supplementary Table 4 ... header text ..."])
    sh.append(["FSTF RNAi", "Cell-type markers tested", "", "", ""])
    sh.append(["dd22331", "dd29413", "dd3524", "dd210", "dd1248"])
    sh.append(["dd11150", "dd15610", "dd8060", "dd28465", None])


def _write_corr(sh: openpyxl.worksheet.worksheet.Worksheet) -> None:
    # Row 1/2/3 are title text, row 4 is the real header, data row 5+
    for _ in range(3):
        sh.append([" Supplementary Table 5 ... title text ..."])
    sh.append(["TF1", "TF2", "X1 Correlation", "G0 Correlation", "G0 Cluster"])
    sh.append(["otp", "elf", 0.06, 0.41, 1])
    sh.append(["otp", "tcf-1", 0.0, 0.17, 1])


def _make_mmc4(xlsx_path: Path) -> None:
    wb = _make_empty_workbook()
    sh = wb.active
    sh.title = "TF"
    _write_tf_catalog(sh)
    wb.save(xlsx_path)


def _make_mmc5(xlsx_path: Path) -> None:
    wb = _make_empty_workbook()
    sh = wb.active
    sh.title = "Sheet1"
    _write_rnai(sh)
    wb.save(xlsx_path)


def _make_mmc6(xlsx_path: Path) -> None:
    wb = _make_empty_workbook()
    sh = wb.active
    sh.title = "Sheet1"
    _write_corr(sh)
    wb.save(xlsx_path)


# ---------------------------------------------------------------------------
# Kings readers
# ---------------------------------------------------------------------------


def test_read_tf_catalog_parses_boolean_flags(tmp_path: Path) -> None:
    p = tmp_path / "mmc4.xlsx"
    _make_mmc4(p)
    df = king.read_tf_catalog(p)
    assert list(df.columns) == [
        "gene_id", "human_hit_acc", "human_hit", "e_value",
        "planarian_name", "percent_identity", "is_tf", "is_prior_fstf",
    ]
    assert df["is_tf"].tolist() == [True, False]
    assert df["is_prior_fstf"].tolist() == [True, False]
    assert df["gene_id"].iloc[0] == "dd_Smed_v6_21801_0_1"


def test_read_rnai_table_returns_long_format(tmp_path: Path) -> None:
    p = tmp_path / "mmc5.xlsx"
    _make_mmc5(p)
    df = king.read_rnai_table(p)
    # First xlsx row: dd22331 → 4 markers; second row: dd11150 → 3 markers (one cell None)
    assert set(df.columns) == {"fstf_rnai", "marker"}
    assert "dd22331" in set(df["fstf_rnai"])
    assert "dd1248" in set(df["marker"])
    assert (df["fstf_rnai"] == "dd22331").sum() == 4
    assert (df["fstf_rnai"] == "dd11150").sum() == 3


def test_read_correlations_round_trips(tmp_path: Path) -> None:
    p = tmp_path / "mmc6.xlsx"
    _make_mmc6(p)
    df = king.read_correlations(p)
    assert list(df.columns) == ["tf1", "tf2", "x1_corr", "g0_corr", "g0_cluster"]
    assert len(df) == 2
    assert df["g0_corr"].iloc[0] == pytest.approx(0.41)


def _add_mmc7_sheet(wb: openpyxl.Workbook, title: str, header_skip: int, header_cols: list[str], subcluster: str) -> None:
    sh = wb.create_sheet(title)
    for _ in range(header_skip):
        sh.append([" title text "])
    sh.append(header_cols)
    sh.append([subcluster, None, None, 1e-115, 5e-62])


def test_read_king_atlas_skips_header_rows(tmp_path: Path) -> None:
    # Build a small mmc7 with all six sheets the reader expects.
    wb = _make_empty_workbook()
    # default empty sheet gets replaced; remove it once we've added real ones.
    wb.remove(wb.active)
    _add_mmc7_sheet(wb, "G0 Progenitor TF Atlas", 4,
                    ["G0 subclusters", "Named Mature Cell Types", "Fincher Clusters", "pou4l-1", "coe"],
                    "neural11_0")
    _add_mmc7_sheet(wb, "G0 Progenitor Pvalues", 3,
                    ["G0 subclusters", "Named Mature Cell Types", "Fincher Clusters", "pv1", "pv2"],
                    "neural11_0")
    _add_mmc7_sheet(wb, "G0 Progenitor Log2FC", 3,
                    ["G0 subclusters", "Named Mature Cell Types", "Fincher Clusters", "v1", "v2"],
                    "neural11_0")
    _add_mmc7_sheet(wb, "X1 TF Atlas", 4,
                    ["X1 subclusters", "Named Cell Types", "f", "x", "y"],
                    "X1_0")
    _add_mmc7_sheet(wb, "X1 Pvalues", 3,
                    ["X1 subclusters", "Named Cell Types", "pv1", "pv2"],
                    "X1_0")
    _add_mmc7_sheet(wb, "X1 Log2FC", 3,
                    ["X1 subclusters", "Named Cell Types", "v1", "v2"],
                    "X1_0")
    wb.save(tmp_path / "mmc7.xlsx")
    out = king.read_king_atlas(tmp_path / "mmc7.xlsx")
    assert "g0_atlas" in out and "g0_pvalues" in out
    assert out["g0_atlas"].columns[0] == "subcluster"
    assert out["g0_atlas"].iloc[0, 0] == "neural11_0"


# ---------------------------------------------------------------------------
# Fincher & Plass readers
# ---------------------------------------------------------------------------


def test_fincher_read_dge_returns_anndata(tmp_path: Path) -> None:
    # Gzip a small DGE-style TSV: genes as rows, cells as columns
    dge = pd.DataFrame(
        {"cell_1": [10, 0, 5], "cell_2": [2, 8, 0], "cell_3": [0, 3, 7]},
        index=["dd_Smed_v4_1", "dd_Smed_v4_2", "dd_Smed_v4_3"],
    )
    p = tmp_path / "dge.txt.gz"
    dge.to_csv(p, sep="\t", compression="gzip")
    adata = fincher.read_dge(p)
    assert isinstance(adata, ad.AnnData)
    assert adata.n_obs == 3 and adata.n_vars == 3
    assert list(adata.var_names) == ["dd_Smed_v4_1", "dd_Smed_v4_2", "dd_Smed_v4_3"]
    assert list(adata.obs_names) == ["cell_1", "cell_2", "cell_3"]


def test_plass_read_h5ad_round_trip(tmp_path: Path) -> None:
    # anndata 0.12 + pandas 3 + pyarrow returns ArrowStringArray for index
    # columns, which anndata can't write without opt-in. Set string storage
    # to plain python for the duration of this test.
    import anndata as _ad
    import pandas as _pd
    _ad.settings.allow_write_nullable_strings = True
    _pd.set_option("mode.string_storage", "python")
    obs = _pd.DataFrame(index=np.array([f"c{i}" for i in range(4)], dtype=object))
    var = _pd.DataFrame(index=np.array([f"g{i}" for i in range(4)], dtype=object))
    adata = ad.AnnData(
        X=np.eye(4, dtype=np.float32),
        obs=obs,
        var=var,
    )
    p = tmp_path / "plass.h5ad"
    adata.write_h5ad(p)
    out = plass.read_plass_matrix(p)
    assert out.n_obs == 4 and out.n_vars == 4


def test_plass_read_txt_gz_matches_fincher_loader(tmp_path: Path) -> None:
    dge = pd.DataFrame(
        {"cell_1": [1, 2], "cell_2": [3, 4]},
        index=["dd_Smed_v6_1", "dd_Smed_v6_2"],
    )
    p = tmp_path / "plass.txt.gz"
    dge.to_csv(p, sep="\t", compression="gzip")
    out = plass.read_plass_matrix(p)
    assert out.n_obs == 2 and out.n_vars == 2
    assert list(out.var_names) == ["dd_Smed_v6_1", "dd_Smed_v6_2"]


def test_plass_read_unsupported_extension_raises(tmp_path: Path) -> None:
    p = tmp_path / "bad.csv"
    p.write_text("not,a,valid,extension")
    with pytest.raises(ValueError):
        plass.read_plass_matrix(p)
