"""NeuralTF pipeline — tiered neural-TF candidate discovery from scRNA-seq atlases.

Usage::
    python -m bioforge.projects.neuraltf.pipeline

Or from the CLI::
    bioforge neuraltf run
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from bioforge.evidence.readers import king as king_reader
from bioforge.evidence import load_bridge
from bioforge.evidence.schema import EvidenceRecord, EvidenceSource
from bioforge.evidence.scoring import EvidenceScorer
from bioforge.evidence.confidence import assign_tiers
from bioforge.evidence.cards import build_cards_for_records, render_cards_markdown


DATA_ROOT = Path.cwd()

_RE_DD_ID = re.compile(r"(dd\D*?\d+)")
_NEURAL_FC_THRESHOLD = 2.0


class NeuralTFPipeline:
    """End-to-end neural TF candidate discovery pipeline.

    Integrates 4 atlases (Fincher, Plass, King, Cui) plus King TF catalog,
    RNAi phenotype table, and neural TF-pair correlations. Produces a
    full ranking and a neural-filtered ranking.

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
        self.bridge_path = self.data_dir / "bridge.csv"
        self.king_atlas_path = self.data_dir / "king_atlas.tsv"
        self.cui_atlas_path = self.data_dir / "cui_atlas_summary.csv"

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
        self.all_records: dict[str, EvidenceRecord] = {}
        self.atlas_membership: dict[str, set[str]] = {}

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
        print("[1/9] Loading datasets...")
        self.adata_fincher = ad.read_h5ad(self.fincher_path)
        print(f"  Fincher: {self.adata_fincher.n_obs} cells x {self.adata_fincher.n_vars} genes (v4)")
        self.adata_plass = ad.read_h5ad(self.plass_path)
        print(f"  Plass:   {self.adata_plass.n_obs} cells x {self.adata_plass.n_vars} genes (v6)")

        if self.subsample:
            for adata, name in [(self.adata_fincher, "Fincher"), (self.adata_plass, "Plass")]:
                if adata.n_obs > self.subsample:
                    sc.pp.subsample(adata, n_obs=self.subsample, random_state=42)
                    print(f"  Subsampled {name} to {adata.n_obs} cells")

    def load_reference_tables(self):
        print("[2/9] Reference tables...")
        self.tf_catalog = pd.read_excel(self.mmc4, sheet_name="TF")
        self.rnai_table = pd.read_excel(self.mmc5, header=None)
        self.correlations = pd.read_excel(self.mmc6, header=None)
        self.tf_ids = set(
            self.tf_catalog.loc[self.tf_catalog["TF?"].notna(), "Gene ID"].astype(str)
        )
        self.tf_ids_norm = self.tf_ids | {tid[:-2] for tid in self.tf_ids if tid.endswith("_1")}
        print(f"  Catalog: {len(self.tf_catalog)} entries ({len(self.tf_ids)} TFs)")
        print(f"  RNAi: {len(self.rnai_table)} rows, Correlations: {len(self.correlations)} pairs")

        # Load Perez 2025 TF classification (MOESM5)
        self.perez_tf_class: dict[str, str] = {}
        perez_path = (
            self.raw_dir / "Supplementary_Data_ Perez_2025"
            / "41467_2025_65712_MOESM5_ESM.xlsx"
        )
        if perez_path.exists():
            try:
                perez = pd.read_excel(perez_path, sheet_name=0, dtype=str, nrows=60000)
                cols = perez.columns.tolist()
                gene_col = cols[0]
                tf_class_col = next((c for c in cols if "TF Class" in c and "Perez" in c), None)
                rbh_col = next((c for c in cols if "1:1" in c and "v6" in c.lower()), None)
                if tf_class_col and rbh_col:
                    for _, r in perez.iterrows():
                        v6 = str(r.get(rbh_col, "")).strip()
                        cls = str(r.get(tf_class_col, "")).strip()
                        if v6 and v6 != "nan" and cls and cls != "nan":
                            self.perez_tf_class[v6] = cls
                    print(f"  Perez TF classification: {len(self.perez_tf_class)} genes")
            except Exception as e:
                print(f"  (Perez TF classification load failed: {e})")

        print("[3/9] Bridge table...")
        self.bridge = load_bridge(self.bridge_path)
        self._enrich_bridge_names()
        print(f"  {len(self.bridge.df)} rows bridged")

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
        """Extract 'dd_Smed_v6_10201_0_1' -> 'dd10201' (or 'dd11150' -> 'dd11150')."""
        if not gene_id or not isinstance(gene_id, str):
            return None
        m = _RE_DD_ID.search(gene_id)
        if not m:
            return None
        raw = m.group(1)
        digits = re.sub(r"\D+", "", raw)
        return f"dd{digits}" if digits else None

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
        print("[4/9] QC + clustering (leiden)...")
        for adata, label in [(self.adata_fincher, "Fincher"), (self.adata_plass, "Plass")]:
            print(f"  {label}: ", end="", flush=True)
            sc.pp.filter_genes(adata, min_cells=3)
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=5000, batch_key=None)
            tf_in_mask = [v for v in adata.var_names if v in self.tf_ids_norm]
            adata.var.loc[tf_in_mask, "highly_variable"] = True
            adata.raw = adata
            hvg = adata[:, adata.var.highly_variable].copy()
            sc.pp.pca(hvg, n_comps=50)
            sc.pp.neighbors(hvg, n_neighbors=10, n_pcs=40)
            sc.tl.leiden(hvg, resolution=0.5)
            adata.obs["leiden"] = hvg.obs["leiden"]
            hvgs = int(adata.var.highly_variable.sum())
            n_cl = adata.obs["leiden"].nunique()
            print(f"cells={adata.n_obs} HVGs={hvgs} leiden={n_cl}")

    # ------------------------------------------------------------------
    # Per-atlas scoring (Fincher, Plass)
    # ------------------------------------------------------------------

    def score_atlases(self):
        print("\n[5/9] Scoring candidates per atlas ...")
        print(f"  {len(self.tf_ids)} TF targets")

        for atlas, atlas_label in [(self.adata_fincher, "fincher"), (self.adata_plass, "plass")]:
            print(f"  {atlas_label}: ", end="", flush=True)
            sc.tl.rank_genes_groups(atlas, "leiden", method="wilcoxon")
            self._score_one_atlas(atlas, atlas_label, atlas.uns["rank_genes_groups"])

        print(f"  Generated records: {len(self.all_records)}")

    def _score_one_atlas(self, adata, atlas_name, result):
        bridge = self.bridge

        if atlas_name == "fincher":
            v6_of = {v: bridge.v4_to_v6(v) for v in adata.var_names}
            score_genes = [v for v, v6 in v6_of.items() if v6 in self.tf_ids_norm]
        else:
            v6_of = {}
            score_genes = [v for v in adata.var_names if v in self.tf_ids_norm]

        n_clusters = len(adata.obs["leiden"].cat.categories)
        clusters = result["names"].dtype.names

        gene_best: dict[str, tuple[float, float]] = {}
        for cl in clusters:
            for g, lfc, pval in zip(result["names"][cl], result["logfoldchanges"][cl], result["pvals"][cl]):
                key = str(g)
                abs_lfc = abs(float(lfc))
                if key not in gene_best or abs_lfc > abs(gene_best[key][0]):
                    gene_best[key] = (abs_lfc, float(pval))

        for gene in score_genes:
            best_l2fc, best_p = gene_best.get(gene, (0.0, 1.0))

            if best_p > 0.05:
                continue

            if atlas_name == "fincher":
                v6_id = bridge.v4_to_v6(gene)
            else:
                v6_id = gene if gene in self.tf_ids else gene + "_1"

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

            # Expression: best (max) per-atlas log2FC across the atlases the
            # gene is enriched in — per-atlas values accumulate (never
            # overwrite an earlier atlas's stronger signal).
            rec.add_score(EvidenceSource.EXPRESSION,
                          max(rec.scores.get(EvidenceSource.EXPRESSION, 0.0),
                              min(1.0, best_l2fc / 5.0)),
                          note=f"log2FC={best_l2fc:.2f},p={best_p:.2g}")

            # specificity [ crude = 1/n_clusters ] — best (most specific)
            # atlas wins; per-atlas values accumulate.
            clusters_here = n_clusters
            rec.add_score(EvidenceSource.SPECIFICITY,
                          max(rec.scores.get(EvidenceSource.SPECIFICITY, 0.0),
                              1.0 / clusters_here if clusters_here > 0 else 0.0),
                          note=f"atlas={atlas_name}")

            self.atlas_membership.setdefault(v6_id, set()).add(atlas_name)

    # ------------------------------------------------------------------
    # King Atlas integration (mmc7, prebuilt TSV)
    # ------------------------------------------------------------------

    def integrate_king_atlas(self):
        print("[6/9] King TF Atlas...")
        if not self.king_atlas_path.exists():
            print("  (missing, skipping)")
            return

        king = pd.read_csv(self.king_atlas_path, sep="\t")

        neural_mask = king["subcluster"].astype(str).str.startswith("neural")
        neural_df = king[neural_mask & (king["log2fc"] >= _NEURAL_FC_THRESHOLD)]

        overall_max = king["log2fc"].max()

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

                feat_expr = min(1.0, fcm / overall_max) if overall_max > 0 else 0.0
                feat_spec = min(1.0, 1.0 / nsub) if nsub > 0 else 0.0

                rec.add_score(EvidenceSource.EXPRESSION,
                              max(rec.scores.get(EvidenceSource.EXPRESSION, 0.0), feat_expr),
                              note=f"king_maxFC={fcm:.2f}")
                rec.add_score(EvidenceSource.SPECIFICITY,
                              max(rec.scores.get(EvidenceSource.SPECIFICITY, 0.0), feat_spec),
                              note=f"king_n_subs={nsub}")

            # Neural-specific signal
            nHits = neural_df[neural_df["v6_id"] == gene_id]
            if len(nHits) > 0:
                rec.add_score(EvidenceSource.NEURAL_ENRICHED, 1.0,
                              note=f"neural_max_l2fc={nHits['log2fc'].max():.2f}")

                unique_ns = nHits["subcluster"].nunique()
                rec.add_score(EvidenceSource.NEURAL_SPECIFICITY,
                              1.0 / unique_ns if unique_ns > 0 else 0.0,
                              note=f"n_neural_subclusters={unique_ns}")

            # King-atlas support only when the gene actually has hits there
            # (a Fincher/Plass-only candidate must not gain king membership).
            if len(hits) > 0:
                self.atlas_membership.setdefault(gene_id, set()).add("king")

    # ------------------------------------------------------------------
    # Cui 2023 Atlas integration (preprocessed CSV)
    # ------------------------------------------------------------------

    def integrate_cui_atlas(self):
        """Integrate Cui 2023 scRNA-seq atlas as 4th evidence source.

        Adds expression and specificity from Cui's 61 cell-type annotations
        across 8 regeneration timepoints. Reuses existing EvidenceSource
        enums (EXPRESSION, SPECIFICITY) — best-atlas-wins semantics.
        """
        print("[6b/9] Cui 2023 Atlas...")
        if not self.cui_atlas_path.exists():
            print(f"  (missing {self.cui_atlas_path}, skipping)")
            return

        cui = pd.read_csv(self.cui_atlas_path)
        print(f"  Loaded {len(cui)} genes from Cui atlas")

        # Cui expression scores are already normalized to [0,1] by the
        # preprocessing script. Neural enrichment/enrichment specificity
        # are also precomputed.
        for _, row in cui.iterrows():
            gene_id = str(row["gene_id"]).strip()
            if not gene_id or gene_id == "nan":
                continue

            if gene_id not in self.all_records:
                self.all_records[gene_id] = EvidenceRecord(gene_id=gene_id)

            rec = self.all_records[gene_id]

            # Expression: best-atlas-wins (same semantics as King)
            cui_expr = float(row.get("expression_score", 0) or 0)
            if cui_expr > 0:
                rec.add_score(
                    EvidenceSource.EXPRESSION,
                    max(rec.scores.get(EvidenceSource.EXPRESSION, 0.0), cui_expr),
                    note=f"cui_fc={row.get('max_fold_change', 0):.2f}",
                )

            # Specificity: best-atlas-wins
            cui_spec = float(row.get("specificity_score", 0) or 0)
            if cui_spec > 0:
                rec.add_score(
                    EvidenceSource.SPECIFICITY,
                    max(rec.scores.get(EvidenceSource.SPECIFICITY, 0.0), cui_spec),
                    note=f"cui_n_types={row.get('n_expressed_types', 0)}",
                )

            # Neural enrichment from Cui (independent of King)
            if row.get("neural_enriched", False):
                # Only upgrade — never downgrade existing neural signal
                existing = rec.scores.get(EvidenceSource.NEURAL_ENRICHED, 0.0)
                if existing < 1.0:
                    rec.add_score(
                        EvidenceSource.NEURAL_ENRICHED, 1.0,
                        note="cui_neural_enriched=True",
                    )

            # Neural specificity from Cui
            cui_nspec = float(row.get("neural_specificity_score", 0) or 0)
            if cui_nspec > 0:
                existing_nspec = rec.scores.get(EvidenceSource.NEURAL_SPECIFICITY, 0.0)
                if cui_nspec > existing_nspec:
                    rec.add_score(
                        EvidenceSource.NEURAL_SPECIFICITY, cui_nspec,
                        note=f"cui_n_neural={row.get('n_neural_expressed', 0)}",
                    )

            self.atlas_membership.setdefault(gene_id, set()).add("cui")

    # ------------------------------------------------------------------
    # RNAi phenotype table (mmc5)
    # ------------------------------------------------------------------

    def integrate_rnai(self):
        print("[8/9] RNAi phenotypes...")
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
            m = _RE_DD_ID.search(val)
            if m:
                targets.add(m.group(1))
        return targets

    # ------------------------------------------------------------------
    # Correlation table (mmc6)
    # ------------------------------------------------------------------

    def integrate_correlations(self):
        print("[9/9] TF pair correlations...")
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
            m = _RE_DD_ID.search(s)
            if m:
                out.add(f"dd{re.sub(r'\\D+', '', m.group(1))}")
            return out

        matched = 0
        for rec in self.all_records.values():
            ids = self._all_ids_for_record(rec)
            if not ids:
                continue
            tf1_ok = data["tf1"].map(normalize)
            tf2_ok = data["tf2"].map(normalize)
            mask = tf1_ok.apply(lambda c: bool(c & ids)) | tf2_ok.apply(lambda c: bool(c & ids))
            sub = data[mask]
            if len(sub) == 0:
                continue
            try:
                x1 = pd.to_numeric(sub["x1_corr"], errors="coerce").mean()
                g0 = pd.to_numeric(sub["g0_corr"], errors="coerce").mean()
            except (ValueError, TypeError):
                continue
            if pd.isna(x1) or pd.isna(g0):
                continue
            gain = max(0.0, float(g0) - float(x1))
            rec.add_score(EvidenceSource.CORRELATION, min(1.0, gain * 3.0),
                          note=f"x1={x1:.2f},g0={g0:.2f}")
            matched += 1
        print(f"  Correlation match: {matched} candidates")

    # ------------------------------------------------------------------
    # Reproducibility (final)
    # ------------------------------------------------------------------

    def assign_reproducibility(self):
        n_atlases = 4  # Fincher, Plass, King, Cui
        for gene_id in self.all_records:
            atlases = self.atlas_membership.get(gene_id, set())
            n = min(len(atlases), n_atlases)
            self.all_records[gene_id].add_score(
                EvidenceSource.REPRODUCIBILITY,
                n / float(n_atlases),
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
                "proof_status": rec.proof_status or "unknown",
                "tier": tier_of.get(rec.gene_id, "low"),
            }
            for src in EvidenceSource:
                row[src.value] = rec.scores.get(src)
            rank_rows.append(row)

        rank_df = pd.DataFrame(rank_rows).sort_values("integrated_score", ascending=False)

        # Neural-specific output
        neural_df = rank_df[
            rank_df["neural_enriched"].notna() | (rank_df["rnai"] > 0)
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
        from typing import Dict
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
    # Main entry
    # ------------------------------------------------------------------

    def run(self):
        self.load_datasets()
        self.load_reference_tables()
        self.run_qc()
        self.score_atlases()
        self.integrate_king_atlas()
        self.integrate_cui_atlas()
        self.integrate_rnai()
        self.integrate_correlations()
        self.assign_reproducibility()
        self.write_outputs()


def main():
    pipe = NeuralTFPipeline()
    pipe.run()


if __name__ == "__main__":
    main()
