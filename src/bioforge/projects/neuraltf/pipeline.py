"""NeuralTF pipeline — tiered neural-TF candidate discovery from scRNA-seq atlases.

Usage::
    python -m bioforge.projects.neuraltf.pipeline

Or from the CLI::
    bioforge neuraltf run
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from statsmodels.stats.multitest import multipletests

from bioforge.evidence import load_bridge
from bioforge.evidence.schema import EvidenceRecord, EvidenceSource
from bioforge.evidence.scoring import EvidenceScorer
from bioforge.evidence.confidence import assign_tiers
from bioforge.evidence.cards import build_cards_for_records, render_cards_markdown


DATA_ROOT = Path(__file__).resolve().parents[4]

# Two gene-ID dialects appear across the King tables and bridge:
#   structured : dd_Smed_v6_11150_0_1   (the numeric gene field is field 3)
#   short      : dd11150                 (King mmc5/mmc6 free-text style)
# The old lazy regex `(dd\D*?\d+)` matched "dd_Smed_v6" on structured IDs
# (digits="6") and silently corrupted short-ID extraction.
_RE_DD_STRUCTURED = re.compile(r"dd_Smed_v[46]_(\d+)")
_RE_DD_SHORT = re.compile(r"\bdd(\d+)\b")
_NEURAL_FC_THRESHOLD = 2.0
_FDR_THRESHOLD = 0.1  # Benjamini-Hochberg q-value threshold
_L2FC_EPS = 1e-9  # pseudocount for true-log2FC means


class NeuralTFPipeline:
    """End-to-end neural TF candidate discovery pipeline.

    Integrates 5 atlases (Fincher 2018, Plass 2018, Cui 2023, King 2024,
    Perez 2025) plus the King TF catalog, RNAi phenotype table,
    neural TF-pair correlations, and Perez 2025 ANANSE regulatory networks.
    Produces a full ranking and a neural-filtered ranking.

    Parameters
    ----------
    data_root : Path
        Root path of the data directory.
    out_dir : Path or None
        Override output directory.
    subsample : int
        Subsample each h5ad to at most this many cells (default 0 =
        keep the complete atlases; a value like 10000 speeds up dev runs).
    """

    def __init__(
        self,
        data_root: Path = DATA_ROOT,
        out_dir: Path | None = None,
        subsample: int = 0,
    ):
        self.data_root = Path(data_root)
        if out_dir is None:
            self.out_dir = self.data_root / "projects" / "NeuralTF" / "runs" / "pipeline_run"
        else:
            self.out_dir = Path(out_dir)
        self.subsample = subsample

        # derived paths
        self.raw_dir = self.data_root / "datasets" / "raw"
        self.proc_dir = self.data_root / "datasets" / "processed"
        self.data_dir = self.data_root / "projects" / "NeuralTF" / "data"

        self.fincher_path = self.proc_dir / "fincher_subsample.h5ad"
        self.plass_path = self.proc_dir / "plass_v6.h5ad"
        self.cui_path = self.proc_dir / "cui_v6.h5ad"
        self.bridge_path = self.data_dir / "bridge.csv"
        self.king_atlas_path = self.data_dir / "king_atlas.tsv"

        king_dir = self.raw_dir / "Supplementary_Data_ King_2024"
        # Try the original Cell Reports (Elsevier) filenames first. If those
        # aren't present, auto-discover mmc4-mmc7.xlsx by suffix (users often
        # rename downloaded supplementary files). The discovery is deterministic
        # — only files matching exactly `mmc4.xlsx`, `mmc5.xlsx` etc. — so
        # accidental co-located files of the same suffix don't collide.
        self.mmc4 = self._resolve_king_xlsx(king_dir, "mmc4")
        self.mmc5 = self._resolve_king_xlsx(king_dir, "mmc5")
        self.mmc6 = self._resolve_king_xlsx(king_dir, "mmc6")
        self.mmc7 = self._resolve_king_xlsx(king_dir, "mmc7")

        # state
        self.tf_catalog: pd.DataFrame | None = None
        self.rnai_table: pd.DataFrame | None = None
        self.correlations: pd.DataFrame | None = None
        self.bridge = None
        self.tf_ids: set[str] = set()
        self.tf_ids_norm: set[str] = set()
        self.tf_ids_king: set[str] = set()  # King mmc4-only seed (HVG forcing)
        self.all_records: dict[str, EvidenceRecord] = {}
        self.atlas_membership: dict[str, set[str]] = {}
        # per-gene best DE (p, lfc) per atlas, for downstream meta-analysis
        self.de_pvals: dict[str, dict[str, tuple[float, float]]] = {}

    @staticmethod
    def _resolve_king_xlsx(king_dir: Path, mmc_name: str) -> Path:
        """Return the path to a King mmcN xlsx by trying the original
        Elsevier filename first, then any `<anything>mmcN.xlsx` in the
        directory, then plain `mmcN.xlsx`. Raises FileNotFoundError only
        when the pipeline actually tries to read the file downstream."""
        # 1) Exact original Elsevier name
        candidate = king_dir / f"1-s2.0-S2211124724001712-{mmc_name}.xlsx"
        if candidate.exists():
            return candidate
        # 2) Glob fallbacks (sorted for deterministic behaviour)
        if king_dir.exists():
            for p in sorted(king_dir.iterdir()):
                # accept "mmc4.xlsx" or "...mmc4.xlsx"
                if p.suffix.lower() == ".xlsx" and (
                    p.stem.lower() == mmc_name
                    or p.stem.lower().endswith(mmc_name)
                ):
                    return p
        # 3) Return the original candidate (default) so the downstream
        # error message references the documented filename.
        return candidate

    # ------------------------------------------------------------------
    # Data Loading
    # ------------------------------------------------------------------

    def load_datasets(self):
        print("[1/10] Loading datasets...")
        self.adata_fincher = ad.read_h5ad(self.fincher_path)
        print(f"  Fincher: {self.adata_fincher.n_obs} cells x {self.adata_fincher.n_vars} genes (v4)")
        self.adata_plass = ad.read_h5ad(self.plass_path)
        print(f"  Plass:   {self.adata_plass.n_obs} cells x {self.adata_plass.n_vars} genes (v6)")
        self.adata_cui = None
        if self.cui_path.exists():
            self.adata_cui = ad.read_h5ad(self.cui_path)
            print(f"  Cui:     {self.adata_cui.n_obs} cells x {self.adata_cui.n_vars} genes (v6)")
        else:
            print(f"  Cui:     (missing {self.cui_path}, skipping)")

        if self.subsample:
            for adata, name in [(self.adata_fincher, "Fincher"), (self.adata_plass, "Plass")]:
                if adata.n_obs > self.subsample:
                    sc.pp.subsample(adata, n_obs=self.subsample, random_state=42)
                    print(f"  Subsampled {name} to {adata.n_obs} cells")

    def load_reference_tables(self):
        print("[2/10] Reference tables...")
        self.tf_catalog = pd.read_excel(self.mmc4, sheet_name="TF")
        self.rnai_table = pd.read_excel(self.mmc5, header=None)
        self.correlations = pd.read_excel(self.mmc6, header=None)
        self.tf_ids = set(
            self.tf_catalog.loc[self.tf_catalog["TF?"].notna(), "Gene ID"].astype(str)
        )
        # King-only seed for HVG forcing: the 14k master catalog would
        # override the data-driven 5,000-HVG selection almost entirely and
        # cluster the atlases in TF space (biasing DE + BH-FDR geometry).
        # The master catalog still seeds records (below) but never HVGs.
        self.tf_ids_king = set(self.tf_ids)
        self.tf_ids_norm = self.tf_ids | {tid[:-2] for tid in self.tf_ids if tid.endswith("_1")}
        print(f"  King catalog: {len(self.tf_catalog)} entries ({len(self.tf_ids)} TFs)")

        # Expand candidate seed with the unified master TF catalog (King + Perez MOESM5).
        # Without this, any TF annotated only in Perez 2025 (up to ~14k genes) would
        # never be tested for de novo cluster DE in Fincher, Plass, or Cui.
        # NOTE: the expanded set seeds records only — run_qc forces just the
        # King mmc4 subset into HVGs (self.tf_ids_king).
        master_path = self.data_dir / "master_tf_catalog.csv"
        if master_path.exists():
            try:
                master = pd.read_csv(master_path, dtype=str)
                if "v6_id" in master.columns:
                    perez_ids = set(master["v6_id"].dropna().str.strip()) - {"", "nan"}
                    n_before = len(self.tf_ids)
                    self.tf_ids |= perez_ids
                    self.tf_ids_norm = (
                        self.tf_ids
                        | {tid[:-2] for tid in self.tf_ids if tid.endswith("_1")}
                    )
                    print(
                        f"  Expanded TF seed: {n_before} (King) -> {len(self.tf_ids)} IDs "
                        f"(+{len(self.tf_ids) - n_before} from Perez MOESM5 / master catalog)"
                    )
            except Exception as e:
                print(f"  (master_tf_catalog expansion failed: {e}; using King seed only)")
        else:
            print(f"  (master_tf_catalog not found at {master_path}; using King seed only)")
        print(f"  RNAi: {len(self.rnai_table)} rows, Correlations: {len(self.correlations)} pairs")

        # Load Perez 2025 TF classification (preprocessed CSV)
        self.perez_tf_class: dict[str, str] = {}
        perez_csv = self.data_dir / "perez_tf_summary.csv"
        if perez_csv.exists():
            try:
                perez = pd.read_csv(perez_csv, dtype=str)
                pcols = perez.columns.tolist()
                if "v6_id" not in pcols or "tf_class" not in pcols:
                    print(f"  WARNING: perez_tf_summary.csv missing expected columns "
                          f"(have {pcols[:5]}..., need 'v6_id' and 'tf_class')")
                else:
                    for _, r in perez.iterrows():
                        v6 = str(r.get("v6_id", "")).strip()
                        cls = str(r.get("tf_class", "")).strip()
                        if v6 and v6 != "nan" and self._valid_perez_class(cls):
                            self.perez_tf_class[v6] = cls
                    print(f"  Perez TF classification: {len(self.perez_tf_class)} genes (from preprocessed CSV)")
            except Exception as e:
                print(f"  (Perez TF classification load failed: {e})")
        else:
            # Fallback: try loading raw MOESM5
            perez_path = (
                self.raw_dir / "Supplementary_Data_ Perez_2025"
                / "41467_2025_65712_MOESM5_ESM.xlsx"
            )
            if perez_path.exists():
                try:
                    perez = pd.read_excel(perez_path, sheet_name=0, dtype=str, nrows=60000)
                    cols = perez.columns.tolist()
                    tf_class_col = next((c for c in cols if "TF Class" in c and "Perez" in c), None)
                    rbh_col = next((c for c in cols if "1:1" in c and "v6" in c.lower()), None)
                    if not tf_class_col:
                        print(f"  WARNING: MOESM5 missing 'TF Class' column. Available: {cols[:8]}...")
                    if not rbh_col:
                        print(f"  WARNING: MOESM5 missing '1:1 v6' column. Available: {cols[:8]}...")
                    if tf_class_col and rbh_col:
                        for _, r in perez.iterrows():
                            v6 = str(r.get(rbh_col, "")).strip()
                            cls = str(r.get(tf_class_col, "")).strip()
                            if v6 and v6 != "nan" and self._valid_perez_class(cls):
                                self.perez_tf_class[v6] = cls
                        print(f"  Perez TF classification: {len(self.perez_tf_class)} genes (from raw MOESM5)")
                except Exception as e:
                    print(f"  (Perez TF classification load failed: {e})")

        print("[3/10] Bridge table...")
        self.bridge = load_bridge(self.bridge_path)
        self._enrich_bridge_names()
        print(f"  {len(self.bridge.df)} rows bridged")

        # Mapping quality report
        try:
            from bioforge.projects.neuraltf.smapping import mapping_stats
            stats = mapping_stats()
            rs = stats["rosetta"]
            m5 = stats["moesm5"]
            print("  Mapping QC:")
            print(f"    Rosetta: {rs['total_smed']} SMED -> {rs['total_v6']} v6 "
                  f"(1-to-many: SMED={rs['smed_one_to_many']}, v6={rs['v6_one_to_many']})")
            print(f"    MOESM5: {m5['total_h1smcg']} h1SMcG -> {m5['total_v6_all']} v6 "
                  f"(h1SMcG->v6 rate: {m5['rate_h1smcg_to_v6']:.2%}, "
                  f"v6->h1SMcG rate: {m5['rate_v6_to_h1smcg']:.2%})")
        except Exception as e:
            print(f"  (Mapping QC skipped: {e})")

    def _enrich_bridge_names(self):
        """Backfill empty bridge.gene_name from mmc4 GenBank names.

        The bridge.csv often has an empty gene_name column even though
        mmc4 holds Planarian GenBank Gene Name for each Gene ID.
        We build a lookup from mmc4 and fill blanks.
        """
        cat = self.tf_catalog
        if cat is None:
            return
        name_col = None
        for c in cat.columns:
            c_low = c.lower()
            if "genbank" in c_low or "gene bank" in c_low or "planarian" in c_low:
                name_col = c
                break
        if name_col is None:
            return
        id_col = None
        for c in cat.columns:
            if c.strip().lower() == "gene id":
                id_col = c
                break
        if id_col is None:
            return

        genbank_of = {}
        for _, r in cat.iterrows():
            gid = str(r[id_col]).strip()
            gname = str(r[name_col]).strip()
            if gname and gname.lower() != "nan" and len(gname) < 50:
                genbank_of[gid] = gname

        df = self.bridge.df
        mask = df["gene_name"].isna() | (df["gene_name"].astype(str).str.strip() == "")
        for idx in df[mask].index:
            v6 = str(df.at[idx, "v6_id"]).strip()
            if v6 in genbank_of:
                df.at[idx, "gene_name"] = genbank_of[v6]

    # ------------------------------------------------------------------
    # Name resolution helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _short_id(gene_id: str) -> str | None:
        """Extract 'dd_Smed_v6_10201_0_1' or 'dd11150' -> 'dd11150'.

        Tries the structured dd_Smed_v[46]_<digits> field first; falls back
        to a short dd<digits> token. Returns None when neither pattern is
        present (non-dd identifiers are never coerced).
        """
        if not gene_id or not isinstance(gene_id, str):
            return None
        m = _RE_DD_STRUCTURED.search(gene_id)
        if m:
            return f"dd{m.group(1)}"
        m = _RE_DD_SHORT.search(gene_id)
        if m:
            return f"dd{m.group(1)}"
        return None

    def _all_ids_for_record(self, record: EvidenceRecord) -> set[str]:
        """Collect every string someone might use to name this gene."""
        ids: set[str] = set()
        if record.gene_name:
            ids.add(record.gene_name)
        sid = self._short_id(record.gene_id)
        if sid:
            ids.add(sid)
        if self.bridge:
            v4 = self.bridge.v6_to_v4(record.gene_id)
            if v4:
                v4_sid = self._short_id(v4)
                if v4_sid:
                    ids.add(v4_sid)
        return ids

    # ------------------------------------------------------------------
    # QC & clustering
    # ------------------------------------------------------------------

    def run_qc(self):
        print("[4/10] QC + clustering (leiden)...")
        atlases = [(self.adata_fincher, "Fincher"), (self.adata_plass, "Plass")]
        if self.adata_cui is not None:
            atlases.append((self.adata_cui, "Cui"))
        for adata, label in atlases:
            print(f"  {label}: ", end="", flush=True)
            sc.pp.filter_cells(adata, min_counts=1)
            sc.pp.filter_genes(adata, min_cells=3)
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            print("norm+log ", end="", flush=True)
            sc.pp.highly_variable_genes(adata, n_top_genes=5000, batch_key=None)
            # Force only the King mmc4 TF catalog (418 genes) into the HVG
            # panel so known TFs enter the clustering; the 14k master catalog
            # would swamp the data-driven HVG selection and cluster in pure
            # TF space (biased DE null + BH-FDR geometry).
            tf_in_mask = [v for v in adata.var_names if v in self.tf_ids_king]
            adata.var.loc[tf_in_mask, "highly_variable"] = True
            adata.raw = adata
            hvg = adata[:, adata.var.highly_variable].copy()
            sc.pp.pca(hvg, n_comps=50)
            sc.pp.neighbors(hvg, n_neighbors=10, n_pcs=40)
            print("neighbors ", end="", flush=True)
            sc.tl.leiden(hvg, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
            adata.obs["leiden"] = hvg.obs["leiden"]
            hvgs = int(adata.var.highly_variable.sum())
            n_cl = adata.obs["leiden"].nunique()
            print(f"leiden={n_cl}")

    # ------------------------------------------------------------------
    # Per-atlas scoring (Fincher, Plass, Cui)
    # ------------------------------------------------------------------

    def score_atlases(self):
        print("\n[5/10] Scoring candidates per atlas ...")
        print(f"  {len(self.tf_ids)} TF targets")

        atlases = [(self.adata_fincher, "fincher"), (self.adata_plass, "plass")]
        if self.adata_cui is not None:
            atlases.append((self.adata_cui, "cui"))
        for atlas, atlas_label in atlases:
            print(f"  {atlas_label}: wilcoxon DE ", end="", flush=True)
            t0 = time.time()
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=pd.errors.PerformanceWarning)
                warnings.simplefilter("ignore", category=FutureWarning)
                sc.tl.rank_genes_groups(atlas, "leiden", method="wilcoxon")
            print(f"({time.time()-t0:.0f}s) scoring ", end="", flush=True)
            t0 = time.time()
            self._score_one_atlas(atlas, atlas_label, atlas.uns["rank_genes_groups"])
            print(f"({time.time()-t0:.0f}s)")


        print(f"  Generated records: {len(self.all_records)}")

    # ------------------------------------------------------------------
    # True log2FC helper
    # ------------------------------------------------------------------

    @staticmethod
    def _cluster_log2fc(adata, gene: str, cluster: str | None) -> float:
        """True log2 fold-change of ``gene`` in ``cluster`` vs all other cells.

        Computed from linear-space normalized expression (``adata.raw``
        holds log1p data; we undo it with expm1, restoring the normalized
        counts-per-10k scale) so the value is a real log2FC, directly
        comparable with King's mmc7 log2FC. A small pseudocount keeps the
        ratio finite for zero-mean denominators.
        """
        if cluster is None:
            return 0.0
        try:
            gene_idx = adata.var_names.get_loc(gene)
        except (KeyError, ValueError):
            return 0.0
        labels = adata.obs["leiden"].astype(str).values
        src = adata.raw if adata.raw is not None else adata
        col = src.X[:, gene_idx]
        vals = np.asarray(col.todense()).ravel() if hasattr(col, "todense") \
            else np.asarray(col).ravel()
        # undo log1p -> linear normalized counts
        lin = np.expm1(vals.astype(np.float64))
        lin = np.where(lin < 0, 0.0, lin)
        in_cl = labels == str(cluster)
        n_in, n_out = int(in_cl.sum()), int((~in_cl).sum())
        if n_in == 0 or n_out == 0:
            return 0.0
        mean_in = float(lin[in_cl].mean())
        mean_out = float(lin[~in_cl].mean())
        fc = (mean_in + _L2FC_EPS) / (mean_out + _L2FC_EPS)
        return float(np.log2(fc)) if fc > 0 else 0.0

    def _score_one_atlas(self, adata, atlas_name, result):
        bridge = self.bridge

        if atlas_name == "fincher":
            v6_of = {v: bridge.v4_to_v6(v) for v in adata.var_names}
            score_genes = [v for v, v6 in v6_of.items() if v6 in self.tf_ids_norm]
        else:
            v6_of = {}
            score_genes = [v for v in adata.var_names if v in self.tf_ids_norm]

        clusters = result["names"].dtype.names

        # Collect all p-values per gene per cluster for global BH-FDR correction
        all_pvals = []
        gene_cluster_keys = []
        for cl in clusters:
            for g, lfc, pval in zip(
                result["names"][cl],
                result["logfoldchanges"][cl],
                result["pvals"][cl],
            ):
                all_pvals.append(float(pval))
                gene_cluster_keys.append((str(g), cl))

        # Benjamini-Hochberg FDR correction (global, across all genes × clusters)
        if all_pvals:
            _, qvals, _, _ = multipletests(all_pvals, alpha=_FDR_THRESHOLD, method="fdr_bh")
        else:
            qvals = []

        # Pre-build per-cluster index lookups to avoid O(N²) scans
        cluster_lookups: dict[str, dict[str, int]] = {}
        for cl in clusters:
            cluster_lookups[cl] = {str(g): i for i, g in enumerate(result["names"][cl])}

        # gene_best: (pos_lfc, pval, qval) — one-tailed: only upregulated signal counts.
        # A repressed gene (lfc < 0) is set to 0.0 so it never inflates the
        # expression score. The best cluster is chosen among SIGNIFICANTLY
        # upregulated clusters (q ≤ threshold); a gene whose largest lfc is
        # non-significant but which has a significant cluster elsewhere is
        # scored on the significant evidence, not dropped entirely.
        gene_best: dict[str, tuple[float, float, float]] = {}
        # gene_best_cluster: argmax significant cluster per gene (for true log2FC)
        gene_best_cluster: dict[str, str] = {}
        # gene_sig_clusters: clusters where the gene is significantly upregulated
        # (qval ≤ threshold AND pos_lfc > 0). Used to compute per-gene specificity.
        gene_sig_clusters: dict[str, set[str]] = {}

        for (key, cl), qval in zip(gene_cluster_keys, qvals):
            idx = cluster_lookups[cl].get(key)
            if idx is None:
                continue
            # One-tailed: only positive fold-changes count as activation evidence
            pos_lfc = max(0.0, float(result["logfoldchanges"][cl][idx]))
            pval = float(result["pvals"][cl][idx])
            if qval <= _FDR_THRESHOLD and pos_lfc > 0:
                # Track all significantly upregulated clusters for specificity
                gene_sig_clusters.setdefault(key, set()).add(cl)
                # Track best significant cluster per gene (max lfc among hits)
                if key not in gene_best or pos_lfc > gene_best[key][0]:
                    gene_best[key] = (pos_lfc, pval, qval)
                    gene_best_cluster[key] = cl
                # per-atlas DE p-value for downstream meta-analysis
                prev = self.de_pvals.get(key, {}).get(atlas_name, (1.0, 0.0))
                if pval < prev[0]:
                    self.de_pvals.setdefault(key, {})[atlas_name] = (pval, pos_lfc)

        for gene in score_genes:
            best_l2fc, best_p, best_q = gene_best.get(gene, (0.0, 1.0, 1.0))

            if best_q > _FDR_THRESHOLD or best_l2fc <= 0:
                continue

            if atlas_name == "fincher":
                v6_id = bridge.v4_to_v6(gene)
            else:
                if gene in self.tf_ids:
                    v6_id = gene
                elif gene + "_1" in self.tf_ids:
                    v6_id = gene + "_1"
                else:
                    continue  # skip genes not in TF catalog

            gn_from_bridge = bridge.v6_to_name(v6_id) if v6_id else None

            if not v6_id:
                continue

            if v6_id not in self.all_records:
                self.all_records[v6_id] = EvidenceRecord(
                    gene_id=v6_id,
                    gene_name=gn_from_bridge,
                )
            rec = self.all_records[v6_id]
            if gn_from_bridge and not rec.gene_name:
                rec.gene_name = gn_from_bridge

            # True log2FC from linear-space cluster means (Scanpy's
            # `logfoldchanges` are a difference of log1p values, not a real
            # log2 fold change). Computed on the best significant cluster
            # so it is directly comparable with King's log2FC.
            true_l2fc = self._cluster_log2fc(adata, gene, gene_best_cluster.get(gene))

            # Expression: best (max) positive log2FC across atlases.
            # Uses one-tailed scoring: min(1, pos_lfc / 5.0) so that a 32-fold
            # upregulation (log2FC=5) saturates the score at 1.0. Repressed
            # genes (lfc<0) get 0 and do not enter this branch.
            rec.add_score(
                EvidenceSource.EXPRESSION,
                max(rec.scores.get(EvidenceSource.EXPRESSION, 0.0),
                    min(1.0, true_l2fc / 5.0)),
                note=f"log2FC={true_l2fc:.2f},p={best_p:.2g},q={best_q:.2g}",
            )

            # Specificity: 1 / n_sig_clusters (clusters where gene is significantly
            # upregulated). A TF expressed in 1 of 16 Leiden clusters scores 1.0;
            # one expressed in all 16 scores 0.0625. This correctly rewards narrow,
            # cell-type-specific activation rather than broad housekeeping expression.
            n_sig = len(gene_sig_clusters.get(gene, set()))
            rec.add_score(
                EvidenceSource.SPECIFICITY,
                max(rec.scores.get(EvidenceSource.SPECIFICITY, 0.0),
                    1.0 / n_sig if n_sig > 0 else 0.0),
                note=f"atlas={atlas_name},n_sig_clusters={n_sig}",
            )

            self.atlas_membership.setdefault(v6_id, set()).add(atlas_name)

    # ------------------------------------------------------------------
    # King Atlas integration (mmc7, prebuilt TSV)
    # ------------------------------------------------------------------

    def integrate_king_atlas(self):
        print("[6/10] King TF Atlas...")
        if not self.king_atlas_path.exists():
            print("  (missing, skipping)")
            return

        king = pd.read_csv(self.king_atlas_path, sep="\t")

        neural_mask = king["subcluster"].astype(str).str.startswith("neural")
        neural_df = king[neural_mask & (king["log2fc"] >= _NEURAL_FC_THRESHOLD)]

        # Expression normalization uses the same 5.0 divisor as the atlas DE
        # path (Fincher/Plass/Cui). log2FC=5.0 = 32-fold upregulation, which
        # maps the physiological TF enrichment range to [0, 1].
        EXPR_CAP = 5.0

        # Pre-compute total subcluster counts for fractional breadth specificity.
        # Fractional breadth: spec = 1 - (n-1)/(N-1) gives 1.0 for n=1 (maximally
        # specific), smoothly declining to 0.0 for n=N (expressed everywhere).
        # This avoids the steep 50% penalty cliff of 1/n when going from n=1→2.
        # Units: (compartment, subcluster) PAIRS on both sides of the fraction —
        # 5 subcluster names (e.g. intestine_*) exist in both G0 and X1, so
        # counting bare names in the denominator against pairs in the numerator
        # mixed units (175 names vs 180 pairs).
        N_KING_TOTAL = max(king.groupby(["compartment", "subcluster"]).ngroups, 2)
        N_NEURAL_TOTAL = max(
            king[neural_mask].groupby(["compartment", "subcluster"]).ngroups, 2
        )

        # Seed records for all TFs enriched in neural G0 compartments.
        # Fincher/Plass subsampling may miss neuron-specific TFs (known
        # from King 2024: neural fate is established post-mitotically and
        # those TFs are poorly resolved in whole-animal Drop-seq at 10K
        # cells). King's G0 progenitor data provides the gold standard.
        # king gene_name is a clean symbol (e.g. "otp", "pax6B") whereas the
        # bridge name is the long mmc4 GenBank description; prefer King's.
        king_name_of = {}
        for _, nr in king.iterrows():
            g = str(nr["gene_name"]).strip()
            if g and g.lower() != "nan" and "transcription factor" not in g.lower():
                king_name_of[str(nr["v6_id"])] = g

        for _, nr in neural_df.iterrows():
            v6_id = nr["v6_id"]
            if v6_id not in self.all_records:
                gn = king_name_of.get(v6_id) or self.bridge.v6_to_name(v6_id)
                self.all_records[v6_id] = EvidenceRecord(
                    gene_id=v6_id, gene_name=gn,
                )
            self.atlas_membership.setdefault(v6_id, set()).add("king")

        for gene_id, rec in self.all_records.items():
            if gene_id in king_name_of and not rec.gene_name:
                rec.gene_name = king_name_of[gene_id]
            # Gene-level push (all of King atlas)
            hits = king[king["v6_id"] == gene_id]
            if len(hits) > 0:
                fcm = hits["log2fc"].max()
                nsub = hits.groupby(["compartment", "subcluster"]).ngroups

                feat_expr = min(1.0, fcm / EXPR_CAP)
                # Fractional breadth specificity: 1 - (n-1)/(N-1)
                # Smooth and biologically graded; avoids the steep 50% cliff of 1/n.
                feat_spec = 1.0 - (nsub - 1) / (N_KING_TOTAL - 1) if nsub > 0 else 0.0
                feat_spec = max(0.0, min(1.0, feat_spec))

                rec.add_score(EvidenceSource.EXPRESSION,
                              max(rec.scores.get(EvidenceSource.EXPRESSION, 0.0), feat_expr),
                              note=f"king_l2fc={fcm:.2f}")
                rec.add_score(EvidenceSource.SPECIFICITY,
                              max(rec.scores.get(EvidenceSource.SPECIFICITY, 0.0), feat_spec),
                              note=f"king_n_subs={nsub},N_total={N_KING_TOTAL}")

            # Neural-specific signal
            nHits = neural_df[neural_df["v6_id"] == gene_id]
            if len(nHits) > 0:
                rec.add_score(EvidenceSource.NEURAL_ENRICHED, 1.0,
                              note=f"neural_max_l2fc={nHits['log2fc'].max():.2f}")

                unique_ns = nHits.groupby(["compartment", "subcluster"]).ngroups
                # Fractional breadth for neural subclusters specifically
                # (same pair units as N_NEURAL_TOTAL).
                neural_spec = 1.0 - (unique_ns - 1) / (N_NEURAL_TOTAL - 1) if unique_ns > 0 else 0.0
                neural_spec = max(0.0, min(1.0, neural_spec))
                rec.add_score(EvidenceSource.NEURAL_SPECIFICITY,
                              neural_spec,
                              note=f"n_neural_subclusters={unique_ns},N_neural={N_NEURAL_TOTAL}")

            # King-atlas support only when the gene actually has hits there
            # (a Fincher/Plass-only candidate must not gain king membership).
            if len(hits) > 0:
                self.atlas_membership.setdefault(gene_id, set()).add("king")

    # ------------------------------------------------------------------
    # Perez 2025 TF lineage classification (EvidenceSource.PEREZ_LINEAGE)
    # ------------------------------------------------------------------

    @staticmethod
    def _valid_perez_class(cls: str) -> bool:
        """True only for a real Perez TF class. MOESM5 uses '-' (and blanks)
        for the 58k non-TF genes — those must never count as TF evidence."""
        return bool(cls) and cls.lower() not in ("nan", "none", "-", "na", "n/a")

    # Neural-relevant TF structural families from Perez 2025 classification.
    # These are the EXACT tf_class strings used in perez_tf_summary.csv (case-insensitive
    # substring match). Selected by literature cross-reference: each family includes
    # conserved planarian neural TFs from King 2024 / Plass 2018 annotation.
    #   bHLH          → Atonal, neurogenin, NeuroD (planarian proneural genes)
    #   Homeodomain   → HD-containing: Otp, Emx, Dlx, Prox, NK, IRX, CUT, MEIS, HOX, PBX
    #   LHX           → LIM-homeodomain: Lhx1, Lhx2, Lhx3, Lhx5 (neural commissure TFs)
    #   POU           → Oct-class: brain2, nub-1 (neoblast neural fate TFs)
    #   Forkhead      → Fox: planarian nervous system FoxG, FoxO
    #   PAX           → Pax3/7, Pax6 (neural progenitor specification)
    #   NKX           → Nkx2.1, Nkx6 (ventral neural patterning)
    #   SIX           → Six3, Six6 (anterior neural regionalization)
    #   EGR           → Egr1 (activity-dependent neural gene)
    #   GLI           → Hedgehog signaling, planarian brain patterning
    #   IRX           → Iroquois, planarian brain segment patterning
    #   COE           → Collier/OLF/EBF (planarian neural identity TFs)
    #   HMG           → SOX-class (all SOX factors have HMG DNA-binding domain)
    #   HOX           → Hox cluster, planarian axis/neural identity
    #   ISL           → Islet-class LIM-homeodomain (motor/sensory neuron fate)
    _PEREZ_NEURAL_CLASSES: frozenset[str] = frozenset({
        "bhlh",         # proneural: atonal, neurogenin, neurod
        "homeodomain",  # broad HD family: otp, emx, dlx, prox, nk, meis, hox, pbx, cut
        "lhx",          # LIM-homeodomain: lhx1/2/3/5 neural
        "pou",          # POU-domain: brain2, nub-1
        "forkhead",     # FOX-domain: foxg, foxo neural
        "pax",          # PAX-domain: pax3/7, pax6
        "nkx",          # NK-homeodomain: nkx2.1, nkx6
        "six",          # SIX-domain: six3, six6 (anterior neural)
        "egr",          # EGR zinc finger (activity-dependent neural)
        "gli",          # GLI hedgehog effectors (brain patterning)
        "irx",          # Iroquois homeodomain (neural segment identity)
        "coe",          # Collier/EBF/OLF (neural identity TFs)
        "hmg",          # HMG-box: SOX factors (neuronal differentiation)
        "hox",          # HOX cluster (axial/neural identity)
        "isl",          # Islet LIM-homeodomain (motor/sensory neuron fate)
    })

    def integrate_perez(self) -> None:
        """Score all records using Perez 2025 TF lineage classification.

        Adds EvidenceSource.PEREZ_LINEAGE to every record:
          - 1.0  : TF class maps to a known neural-relevant family (bHLH, Homeobox, POU, etc.)
          - 0.5  : Gene has a TF class in Perez but not neural-specific
          - 0.0  : Gene absent from Perez MOESM5 or no TF class recorded
        """
        print("[7/10] Perez TF classification...")
        if not self.perez_tf_class:
            print("  (Perez TF classification empty, skipping)")
            return

        matched_neural = 0
        matched_other = 0
        unmatched = 0

        for v6_id, rec in self.all_records.items():
            cls = self.perez_tf_class.get(v6_id, "")
            if not cls and v6_id.endswith("_1"):
                # Also try without the _1 suffix variant (strip the exact
                # suffix; rstrip("_1") would strip any trailing {_,1} chars)
                cls = self.perez_tf_class.get(v6_id.removesuffix("_1"), "")
            if not cls:
                rec.add_score(EvidenceSource.PEREZ_LINEAGE, 0.0,
                              note="absent_from_perez")
                unmatched += 1
                continue
            cls_lower = cls.lower()
            is_neural = any(nf in cls_lower for nf in self._PEREZ_NEURAL_CLASSES)
            score = 1.0 if is_neural else 0.5
            rec.add_score(EvidenceSource.PEREZ_LINEAGE, score,
                          note=f"perez_class={cls}")
            # Gate atlas membership: only genes with a confirmed TF class (score ≥ 0.5)
            # count as supported by Perez atlas. Absent genes (score = 0.0) would
            # otherwise inflate the reproducibility denominator incorrectly.
            if score >= 0.5:
                self.atlas_membership.setdefault(v6_id, set()).add("perez")
            if is_neural:
                matched_neural += 1
            else:
                matched_other += 1

        print(f"  Perez scores: {matched_neural} neural-class, "
              f"{matched_other} other-class, {unmatched} absent")

    # ------------------------------------------------------------------
    # Perez 2025 ANANSE regulatory influence (MOESM19)
    # ------------------------------------------------------------------

    def integrate_perez_influence(self) -> None:
        """Score all records using Perez 2025 ANANSE regulatory influence.

        Adds EvidenceSource.PEREZ_INFLUENCE to every record found as a
        regulatory factor in the neuron fate ANANSE network (MOESM19).
        The influence_score is a normalised 0-1 rank per fate, where 1.0
        means the TF has the highest regulatory influence in that fate.
        """
        print("[8/10] Perez 2025 ANANSE influence...")
        moesm19 = (
            self.raw_dir / "Supplementary_Data_ Perez_2025"
            / "41467_2025_65712_MOESM19_ESM.xlsx"
        )
        if not moesm19.exists():
            print(f"  (MOESM19 not found at {moesm19}, skipping)")
            return

        try:
            xl = pd.ExcelFile(moesm19)
        except Exception as e:
            print(f"  (MOESM19 read failed: {e}, skipping)")
            return

        # Find the neuron fate sheet
        neuron_sheet = None
        for name in xl.sheet_names:
            if "neuron" in name.lower():
                neuron_sheet = name
                break
        if neuron_sheet is None:
            print("  (no neuron fate sheet found in MOESM19, skipping)")
            return

        try:
            infl_df = pd.read_excel(moesm19, sheet_name=neuron_sheet, dtype=str)
        except Exception as e:
            print(f"  (MOESM19 neuron sheet read failed: {e}, skipping)")
            return

        print(f"  MOESM19 neuron sheet: {len(infl_df)} factors")

        # Build h1SMcG -> influence_score lookup
        infl_col = None
        for c in infl_df.columns:
            if "influence_score" in str(c).lower() and "raw" not in str(c).lower():
                infl_col = c
                break
        if infl_col is None:
            print("  (influence_score column not found, skipping)")
            return

        factor_col = infl_df.columns[0]  # 'factor' column (h1SMcG IDs)
        h1_to_influence: dict[str, float] = {}
        for _, row in infl_df.iterrows():
            h1 = str(row[factor_col]).strip()
            try:
                score = float(row[infl_col])
            except (ValueError, TypeError):
                continue
            if h1 and h1 != "nan" and 0.0 <= score <= 1.0:
                h1_to_influence[h1] = score

        print(f"  h1SMcG factors with influence scores: {len(h1_to_influence)}")

        # Map v6 -> h1SMcG and score. RBH-restricted: the collapsed
        # `Similar` mapping claims 14.4k of 25k v6 IDs for >1 h1SMcG, so
        # first-wins picks would attribute influence arbitrarily.
        from bioforge.projects.neuraltf.smapping import batch_v6_to_h1smcg
        v6_ids = list(self.all_records.keys())
        v6_to_h1 = batch_v6_to_h1smcg(v6_ids, rbh_only=True)

        matched = 0
        for v6_id, rec in self.all_records.items():
            h1 = v6_to_h1.get(v6_id)
            if h1 and h1 in h1_to_influence:
                rec.add_score(EvidenceSource.PEREZ_INFLUENCE,
                              h1_to_influence[h1],
                              note=f"neuron_influence={h1_to_influence[h1]:.3f}")
                self.atlas_membership.setdefault(v6_id, set()).add("perez")
                matched += 1

        print(f"  Perez influence match: {matched} candidates")

    # ------------------------------------------------------------------
    # RNAi phenotype table (mmc5)
    # ------------------------------------------------------------------

    def integrate_rnai(self):
        print("[9/10] RNAi phenotypes...")
        if self.rnai_table is None:
            print("  (no table, skipping)")
            return

        rnai_targets = self._build_rna_target_set()
        print(f"  Parsed {len(rnai_targets)} RNAi targets from mmc5")

        matched = 0
        for rec in self.all_records.values():
            ids = self._all_ids_for_record(rec)
            hit = any(x in rnai_targets for x in ids if x)
            rec.add_score(EvidenceSource.RNai, 1.0 if hit else 0.0,
                          note=f"in_mmc5={hit}")
            if hit:
                matched += 1
        print(f"  RNAi match: {matched} candidates")

    def _build_rna_target_set(self) -> set[str]:
        """Parse mmc5 column 0 into a cleaned target set.

        This cleanly handles three formats found in mmc5:
            dd11150
            dd22163 (UNCX)
            fer3l-1
            pax2b
        """
        targets: set[str] = set()
        for _, row in self.rnai_table.iterrows():
            val = str(row.iloc[0]).strip()
            if not val or val == "nan" or val.startswith("Supplementary") or val == "FSTF RNAI":
                continue
            if "transcription factors" in val.lower() or "marker" in val.lower():
                continue
            # Strip parenthetical annotations
            clean = val.split(" (")[0].strip()
            if clean and not clean.startswith("All"):
                targets.add(clean)
            m = _RE_DD_SHORT.search(val) or _RE_DD_STRUCTURED.search(val)
            if m:
                targets.add(f"dd{m.group(1)}")
        return targets

    # ------------------------------------------------------------------
    # Correlation table (mmc6)
    # ------------------------------------------------------------------

    def integrate_correlations(self):
        print("[10/10] TF pair correlations...")
        if self.correlations is None or self.correlations.shape[1] < 4:
            print("  No correlations available")
            return

        data = self.correlations.iloc[4:].copy()
        data.columns = ["tf1", "tf2", "x1_corr", "g0_corr", "g0_cluster"]

        def normalize(val) -> set[str]:
            s = str(val).strip()
            if not s or s.lower() == "nan":
                return set()
            out = {s}
            if " (" in s:
                out.add(s.split(" (", 1)[0].strip())
            m = _RE_DD_STRUCTURED.search(s) or _RE_DD_SHORT.search(s)
            if m:
                out.add(f"dd{m.group(1)}")
            return out

        # Pre-compute normalization once (not per candidate)
        tf1_ok = data["tf1"].map(normalize)
        tf2_ok = data["tf2"].map(normalize)

        matched = 0
        for rec in self.all_records.values():
            ids = self._all_ids_for_record(rec)
            if not ids:
                continue
            mask = tf1_ok.apply(lambda c: bool(c & ids)) | tf2_ok.apply(lambda c: bool(c & ids))
            sub = data[mask]
            if len(sub) == 0:
                continue
            # Joint NaN masking: x1/g0 for the SAME pair row must align.
            # Dropping NaNs independently per column can silently re-pair
            # an x1 value with a different row's g0 value.
            x1_num = pd.to_numeric(sub["x1_corr"], errors="coerce")
            g0_num = pd.to_numeric(sub["g0_corr"], errors="coerce")
            pair = pd.DataFrame({"x1": x1_num, "g0": g0_num}).dropna()
            if pair.empty:
                continue
            # Use max(Δr) across all tested partner pairs rather than mean.
            # A TF may form a tight heterodimeric complex with one specific partner
            # (high Δr) while being uncorrelated with other tested partners (Δr≈0).
            # Averaging dilutes this authentic signal; max captures the best evidence
            # of post-mitotic G0 co-activation for any single functional partnership.
            pair_gains = (pair["g0"] - pair["x1"]).to_numpy()
            best_gain = float(np.max(pair_gains))
            gain = max(0.0, best_gain)
            best_idx = int(np.argmax(pair_gains))
            rec.add_score(
                EvidenceSource.CORRELATION,
                min(1.0, gain * 3.0),
                note=f"best_pair_x1={pair['x1'].iloc[best_idx]:.2f},"
                     f"g0={pair['g0'].iloc[best_idx]:.2f},"
                     f"delta_r={best_gain:.2f},n_pairs={len(pair)}",
            )
            matched += 1
        print(f"  Correlation match: {matched} candidates")

    # ------------------------------------------------------------------
    # Reproducibility (final)
    # ------------------------------------------------------------------

    def assign_reproducibility(self):
        # Reproducibility across the 5 single-cell / TF regulatory atlases:
        # Fincher 2018, Plass 2018, Cui 2023, King 2024, Perez 2025
        n_atlases = 5
        for gene_id in self.all_records:
            atlases = self.atlas_membership.get(gene_id, set())
            self.all_records[gene_id].add_score(
                EvidenceSource.REPRODUCIBILITY,
                len(atlases) / float(n_atlases),
                note=f"atlases={sorted(atlases)}",
            )

    # ------------------------------------------------------------------
    # Outputs
    # ------------------------------------------------------------------

    def write_outputs(self):
        scorer = EvidenceScorer()
        records = list(self.all_records.values())
        tiered = assign_tiers(records, scorer=scorer)
        tiered.sort(key=lambda t: t[2], reverse=True)

        prior_fstf_ids = set(
            self.tf_catalog.loc[
                self.tf_catalog["FSTF?"].astype(str).str.strip().str.lower() == "yes",
                "Gene ID",
            ]
        )

        cards = build_cards_for_records(
            records,
            atlas_membership={
                g: sorted(v) for g, v in self.atlas_membership.items()
            },
            prior_fstf_ids=prior_fstf_ids,
        )
        for card in cards:
            parent = next((r for r in records if r.gene_id == card.gene_id), None)
            if parent:
                card.integrated_score = scorer.integrated_score(parent)

        self.out_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'=' * 70}")
        print(f"{'TIER':<14} {'SCORE':<8} {'GENE':<22} {'STATUS':<22} {'SOURCES'}")
        print(f"{'=' * 70}")
        for rec, tier, score in tiered[:40]:
            st = rec.proof_status or "?"
            srcs = rec.supporting_streams()
            print(f"  {tier.value:<12}  {score:<6.4f}  {(rec.gene_name or rec.gene_id):<20}  {st:<20}  {srcs}")

        # Build a lookup from gene_id -> tier so the CSV carries tier too.
        tier_of = {rec.gene_id: tier.value for rec, tier, _ in tiered}

        rank_rows = []
        for rec in records:
            row = {
                "gene_id": rec.gene_id,
                "gene_name": rec.gene_name or "",
                "integrated_score": scorer.integrated_score(rec),
                "n_streams": rec.supporting_streams(),
                "completeness": round(rec.completeness, 3),
                "proof_status": rec.proof_status or "unknown",
                "tier": tier_of.get(rec.gene_id, "low"),
            }
            for src in EvidenceSource:
                row[src.value] = rec.scores.get(src)
            rank_rows.append(row)

        rank_df = pd.DataFrame(rank_rows).sort_values("integrated_score", ascending=False)

        # Neural-specific output: genes with neural enrichment evidence OR RNAi phenotype.
        # The .fillna(0) guard prevents ValueError when the rnai column is all-NaN
        # (which can happen on small subsampled dev runs where no RNAi targets were scored).
        neural_df = rank_df[
            (rank_df["neural_enriched"].notna() & (rank_df["neural_enriched"] > 0))
            | (rank_df["rnai"].fillna(0) > 0)
        ]
        if len(neural_df) == 0:
            neural_df = rank_df.head(25)

        # Write full
        rpath = self.out_dir / "rank.csv"
        rank_df.to_csv(rpath, index=False)
        print(f"\nFull rank: {rpath} ({len(rank_df)} candidates)")

        # Write neural-filtered
        nr_path = self.out_dir / "rank_neural.csv"
        neural_df.to_csv(nr_path, index=False)
        print(f"Neural rank: {nr_path} ({len(neural_df)} candidates)")

        # Evidence cards
        md = render_cards_markdown(cards)
        cpath = self.out_dir / "evidence_cards.md"
        cpath.write_text(md, encoding="utf-8")
        print(f"Cards:   {cpath} ({len(cards)} cards)")

        # Quick terminal summary
        print("\n=== NEURAL CANDIDATES FOR VALIDATION ===")
        for _, row in neural_df.head(30).iterrows():
            nm = row["gene_name"] or row["gene_id"]
            st = row["proof_status"]
            print(f"  {nm:<20} {row['integrated_score']:.3f}  {st}")

        # JSON
        top_50: list[dict[str, Any]] = []
        for rec, tier, score in tiered[:50]:
            top_50.append({
                "gene_id": rec.gene_id,
                "gene_name": rec.gene_name,
                "score": round(score, 4),
                "tier": tier.value,
                "proof_status": rec.proof_status or "unknown",
                "sources": sorted([k for k in rec.scores]),
            })

        payload = {
            "n_records": len(records),
            "n_cards": len(cards),
            "n_tiered": len(tiered),
            "top_candidates": top_50,
        }
        json_path = self.out_dir / "pipeline_results.json"
        json_path.write_text(json.dumps(payload, indent=2))
        print(f"JSON:    {json_path}")
        print("Done.")

    # ------------------------------------------------------------------
    # Checkpointing helpers
    # ------------------------------------------------------------------

    def _write_checkpoint(self, name: str, df: pd.DataFrame) -> None:
        """Write a pipeline checkpoint file for post-run auditability.

        Parquet is preferred (preserves dtypes). Falls back to CSV if pyarrow
        is not available. Files are written to ``self.out_dir``.
        """
        self.out_dir.mkdir(parents=True, exist_ok=True)
        parquet_path = self.out_dir / f"{name}.parquet"
        csv_path = self.out_dir / f"{name}.csv"
        try:
            df.to_parquet(parquet_path, index=False)
            print(f"  [checkpoint] {parquet_path.name} ({len(df)} rows)")
        except Exception:
            df.to_csv(csv_path, index=False)
            print(f"  [checkpoint] {csv_path.name} ({len(df)} rows)")

    def _checkpoint_atlas_loads(self) -> None:
        rows = []
        for label, adata in [
            ("fincher", self.adata_fincher),
            ("plass", self.adata_plass),
            ("cui", self.adata_cui),
        ]:
            if adata is not None:
                rows.append({
                    "atlas": label,
                    "n_cells": adata.n_obs,
                    "n_genes": adata.n_vars,
                    "has_leiden": "leiden" in adata.obs.columns,
                })
        self._write_checkpoint("checkpoint_01_atlas_loads", pd.DataFrame(rows))

    def _checkpoint_post_qc(self) -> None:
        rows = []
        for label, adata in [
            ("fincher", self.adata_fincher),
            ("plass", self.adata_plass),
            ("cui", self.adata_cui),
        ]:
            if adata is not None:
                n_hvg = int(adata.var.get("highly_variable", pd.Series()).sum()) \
                    if "highly_variable" in adata.var.columns else 0
                rows.append({
                    "atlas": label,
                    "n_cells_post_qc": adata.n_obs,
                    "n_genes_post_qc": adata.n_vars,
                    "n_hvg": n_hvg,
                    "n_leiden_clusters": adata.obs["leiden"].nunique()
                    if "leiden" in adata.obs.columns else 0,
                })
        self._write_checkpoint("checkpoint_02_post_qc", pd.DataFrame(rows))

    def _checkpoint_king_records(self) -> None:
        rows = [{"v6_id": gid, "n_streams": rec.supporting_streams(),
                 "has_neural_enriched": EvidenceSource.NEURAL_ENRICHED in rec.scores}
                for gid, rec in self.all_records.items()]
        self._write_checkpoint("checkpoint_04_king_records", pd.DataFrame(rows))

    def _checkpoint_perez_records(self) -> None:
        rows = [{"v6_id": gid,
                 "perez_score": rec.scores.get(EvidenceSource.PEREZ_LINEAGE, float("nan")),
                 "perez_note": rec.notes.get(EvidenceSource.PEREZ_LINEAGE, "")}
                for gid, rec in self.all_records.items()]
        self._write_checkpoint("checkpoint_05_perez_records", pd.DataFrame(rows))

    def _checkpoint_stream_matrix(self) -> None:
        rows = []
        for gid, rec in self.all_records.items():
            row: dict = {"v6_id": gid, "gene_name": rec.gene_name or ""}
            for src in EvidenceSource:
                row[src.value] = rec.scores.get(src, float("nan"))
            rows.append(row)
        self._write_checkpoint("checkpoint_06_stream_matrix", pd.DataFrame(rows))

    # ------------------------------------------------------------------
    # Main entry
    # ------------------------------------------------------------------

    def _checkpoint_post_scoring(self) -> None:
        """Checkpoint 03: records seeded from atlas DE after score_atlases().

        Also persists per-gene best DE (p, lfc) per atlas for downstream
        meta-analysis (scripts/stats/meta_analytic_pvalue.py).
        """
        rows = [{"v6_id": gid, "n_streams_so_far": rec.supporting_streams(),
                 "has_expression": EvidenceSource.EXPRESSION in rec.scores,
                 "has_specificity": EvidenceSource.SPECIFICITY in rec.scores}
                for gid, rec in self.all_records.items()]
        self._write_checkpoint("checkpoint_03_post_scoring", pd.DataFrame(rows))

        de_rows = []
        for gene, per_atlas in self.de_pvals.items():
            if self.bridge is None:
                v6_id = gene if gene in self.tf_ids else gene + "_1" \
                    if gene + "_1" in self.tf_ids else ""
            else:
                v6_id = self.bridge.v4_to_v6(gene) if gene.startswith("dd_Smed_v4") \
                    else (gene if gene in self.tf_ids else
                          (gene + "_1" if gene + "_1" in self.tf_ids else ""))
            if not v6_id:
                continue
            row = {"v6_id": v6_id, "gene_id_atlas": gene}
            for atlas in ("fincher", "plass", "cui"):
                if atlas in per_atlas:
                    row[f"{atlas}_p"] = per_atlas[atlas][0]
                    row[f"{atlas}_lfc"] = per_atlas[atlas][1]
            de_rows.append(row)
        if de_rows:
            # Aggregate to ONE row per v6_id: different atlases refer to
            # the same v6 gene via different local IDs (v4 vs v6 vs
            # suffix variants), so a naive concat leaves one row per
            # (gene, atlas) and a keep-first dedup would silently drop
            # two of the three atlases' p-values.
            df = pd.DataFrame(de_rows)
            agg = {"gene_id_atlas": "first"}
            for atlas in ("fincher", "plass", "cui"):
                agg[f"{atlas}_p"] = "min"    # best (smallest) p per atlas
                agg[f"{atlas}_lfc"] = "first"
            df = df.groupby("v6_id", as_index=False).agg(agg)
            self._write_checkpoint("de_pvalues", df)

    def run(self):
        self.load_datasets()
        self._checkpoint_atlas_loads()          # checkpoint 01
        self.load_reference_tables()
        self.run_qc()
        self._checkpoint_post_qc()              # checkpoint 02
        self.score_atlases()
        self._checkpoint_post_scoring()         # checkpoint 03
        self.integrate_king_atlas()
        self._checkpoint_king_records()         # checkpoint 04
        self.integrate_perez()
        self._checkpoint_perez_records()        # checkpoint 05
        self.integrate_perez_influence()
        self.integrate_rnai()
        self.integrate_correlations()
        self.assign_reproducibility()
        self._checkpoint_stream_matrix()        # checkpoint 06
        self.write_outputs()


def main():
    pipe = NeuralTFPipeline()
    pipe.run()


if __name__ == "__main__":
    main()
