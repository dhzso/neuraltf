#!/usr/bin/env python
"""ANANSE Full-Candidate GRN Regulatory Network Scan.

Cross-references ALL prioritized TF candidates (from rank.csv) against the
ANANSE-predicted TF-target regulatory network from Perez et al. 2025 (MOESM22).

Unlike `validate_with_perez.py` which only checks the top-10 shortlist, this
script scans every candidate and persists a full edge/target matrix.

Method:
  1. Load all candidates from rank.csv (full set, ~249+).
  2. Load MOESM22 (all fates, not just neuron).
  3. Map candidate v6_ids â†’ h1SMcG via batch_v6_to_h1smcg().
  4. For each candidate found as a TF node in ANANSE:
       - Record fate(s), n_targets, top target symbols, out-degree.
  5. For each candidate found as a target of ANANSE TFs:
       - Record which TFs regulate it and their gene symbols.

Outputs:
  results/ananse_network_full.csv      â€” all candidates with network metrics
  results/ananse_top_regulators.csv    â€” top 20 by out-degree (TF targets)
  results/ananse_network_full.parquet  â€” machine-readable version

Schema of ananse_network_full.csv:
  v6_id, gene_name, h1smcg_id, is_ananse_tf, is_ananse_target,
  fates_as_tf, n_targets_total, n_targets_neuron, top_5_targets,
  regulating_tfs, rank_position, integrated_score, proof_status

Usage:
    python projects/NeuralTF/scripts/ananse_full_scan.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from bioforge.projects.neuraltf.smapping import batch_v6_to_h1smcg   # noqa: E402

RESULTS  = ROOT / "projects" / "NeuralTF" / "results"
RUN      = ROOT / "projects" / "NeuralTF" / "runs" / "pipeline_run"
MOESM22  = (
    ROOT / "datasets" / "raw" / "Supplementary_Data_ Perez_2025"
    / "41467_2025_65712_MOESM22_ESM.xlsx"
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rank_path() -> Path:
    for name in ("rank_all_candidates.csv", "rank.csv"):
        p = RUN / name
        if p.exists():
            return p
    raise FileNotFoundError(
        f"rank.csv or rank_all_candidates.csv not found in {RUN}. "
        "Run the pipeline (scripts/run.py) first."
    )


def load_ananse(path: Path) -> pd.DataFrame:
    """Load all MOESM22 sheets and concatenate."""
    xl = pd.ExcelFile(path)
    sheets = xl.sheet_names
    print(f"  MOESM22 sheets: {sheets}")

    dfs = []
    for sheet in sheets:
        try:
            df = pd.read_excel(path, sheet_name=sheet, dtype=str)
            df["_sheet"] = sheet
            dfs.append(df)
        except Exception as e:
            print(f"  [warn] could not read sheet '{sheet}': {e}")
    if not dfs:
        raise ValueError(f"No readable sheets in {path}")
    combined = pd.concat(dfs, ignore_index=True)
    # Normalise column names
    combined.columns = [str(c).strip() for c in combined.columns]
    return combined


def _col(df: pd.DataFrame, *keywords: str) -> str | None:
    """Find first column whose name contains any of the keywords (case-insensitive)."""
    for kw in keywords:
        for c in df.columns:
            if kw.lower() in c.lower():
                return c
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    if not MOESM22.exists():
        print(f"[WARN] MOESM22 not found: {MOESM22}")
        print("       Download Perez 2025 MOESM22 and place it under:")
        print(f"       {MOESM22.parent}")
        return 1

    rank_path = _rank_path()

    print("== ANANSE Full-Candidate GRN Scan ==")
    print(f"  candidates: {rank_path}")
    print(f"  ANANSE    : {MOESM22.name}")

    # ---------- Load candidates -------------------------------------------
    rank = pd.read_csv(rank_path)
    assert len(rank) > 0, "rank.csv is empty"
    print(f"  Loaded {len(rank)} candidates")

    v6_ids = [str(v) for v in rank["gene_id"] if str(v) != "nan"]
    gene_names = dict(zip(
        rank["gene_id"].astype(str),
        rank.get("gene_name", rank["gene_id"]).fillna("").astype(str),
    ))
    scores = dict(zip(
        rank["gene_id"].astype(str),
        rank.get("integrated_score", pd.Series([0.0] * len(rank))).fillna(0),
    ))
    proof = dict(zip(
        rank["gene_id"].astype(str),
        rank.get("proof_status", pd.Series([""] * len(rank))).fillna(""),
    ))
    rank_pos = {gid: i + 1 for i, gid in enumerate(rank["gene_id"].astype(str))}

    # ---------- Map v6 → h1SMcG -------------------------------------------
    # RBH-restricted: the collapsed Similar column claims 14.4k of 25k v6
    # IDs for >1 h1SMcG, so first-wins picks would attribute edges
    # arbitrarily. The 1:1 RBH column is unambiguous.
    print("  Mapping v6 IDs → h1SMcG (RBH-only) ...")
    v6_to_h1 = batch_v6_to_h1smcg(v6_ids, rbh_only=True)
    h1_to_v6 = {h: v for v, h in v6_to_h1.items() if h}
    n_mapped = sum(1 for h in v6_to_h1.values() if h)
    print(f"  Mapped {n_mapped}/{len(v6_ids)} v6 IDs to h1SMcG (1:1 RBH)")

    # ---------- Load ANANSE -----------------------------------------------
    print(f"  Loading ANANSE network from {MOESM22.name} ...")
    ananse = load_ananse(MOESM22)
    print(f"  ANANSE total interactions: {len(ananse)}")

    tf_col = _col(ananse, "TF (gene ID)", "TF gene ID", "TF_gene_id", "TF")
    tgt_col = _col(ananse, "Target gene (gene ID)", "Target gene", "Target_gene")
    fate_col = _col(ananse, "Fate", "fate", "cell_type")
    tf_sym_col = _col(ananse, "TF (gene symbol)", "TF_symbol")
    tgt_sym_col = _col(ananse, "Target (gene symbol)", "Target_symbol")

    if tf_col is None or tgt_col is None:
        raise ValueError(
            f"Could not identify TF/Target columns in MOESM22. "
            f"Columns present: {list(ananse.columns)}"
        )

    print(f"  TF column: '{tf_col}', Target column: '{tgt_col}'")
    score_col = _col(ananse, "interaction score", "score")
    if fate_col:
        fates = ananse[fate_col].dropna().unique()
        print(f"  Fates: {sorted(fates)}")
    if score_col:
        ananse[score_col] = pd.to_numeric(ananse[score_col], errors="coerce")

    # Per-fate edge totals (for normalization: raw cross-fate edge counts
    # favor fates with more sheets/edges, e.g. phagocytes 2253 vs neuron 1115)
    fate_totals = (
        ananse.dropna(subset=[tf_col])
        .groupby(fate_col)[tf_col].size().to_dict()
        if fate_col else {}
    )

    # Build lookup sets
    ananse_tf_ids: set[str] = set(ananse[tf_col].dropna().astype(str))
    ananse_tgt_ids: set[str] = set(ananse[tgt_col].dropna().astype(str))

    # ---------- Scan each candidate ----------------------------------------
    records = []
    for v6_id in v6_ids:
        h1 = v6_to_h1.get(v6_id, "")

        is_tf = bool(h1 and h1 in ananse_tf_ids)
        is_tgt = bool(h1 and h1 in ananse_tgt_ids)

        # TF metrics
        fates_as_tf: list[str] = []
        n_targets_total = 0
        n_targets_neuron = 0
        top5_targets = ""
        neuron_share = 0.0  # share of a fate's full edge set this TF covers
        if is_tf:
            tf_rows = ananse[ananse[tf_col] == h1]
            n_targets_total = tf_rows[tgt_col].nunique() if tgt_col else len(tf_rows)
            if fate_col:
                fates_as_tf = sorted(tf_rows[fate_col].dropna().unique().tolist())
                neuron_rows = tf_rows[
                    tf_rows[fate_col].str.lower().str.contains("neuron", na=False)
                ]
                n_targets_neuron = (
                    neuron_rows[tgt_col].nunique() if tgt_col else len(neuron_rows)
                )
                # neuron-fate normalized share (per-fate edge-count scale)
                tot = fate_totals.get("neuron", 0)
                if tot > 0:
                    neuron_share = round(len(neuron_rows) / tot, 4)
            else:
                fates_as_tf = ["unknown"]
            # Top targets ordered by interaction score (not file order)
            ordered = tf_rows
            if score_col and score_col in tf_rows.columns:
                ordered = tf_rows.sort_values(score_col, ascending=False)
            sym_src = tgt_sym_col if (tgt_sym_col and tgt_sym_col in ordered.columns) \
                else tgt_col
            if sym_src in ordered.columns:
                seen: list[str] = []
                for x in ordered[sym_src].dropna().astype(str):
                    if x not in seen:
                        seen.append(x)
                    if len(seen) >= 5:
                        break
                top5_targets = "; ".join(seen)

        # Target metrics (which TFs regulate this gene)
        regulating_tfs = ""
        if is_tgt:
            tgt_rows = ananse[ananse[tgt_col] == h1]
            if tf_sym_col and tf_sym_col in tgt_rows.columns:
                reg = tgt_rows[tf_sym_col].dropna().unique()[:5]
            else:
                reg = tgt_rows[tf_col].dropna().unique()[:5]
            regulating_tfs = "; ".join(str(x) for x in reg)

        records.append({
            "v6_id": v6_id,
            "gene_name": gene_names.get(v6_id, ""),
            "h1smcg_id": h1 or "",
            "is_ananse_tf": is_tf,
            "is_ananse_target": is_tgt,
            "fates_as_tf": "|".join(fates_as_tf),
            "n_targets_total": n_targets_total,
            "n_targets_neuron": n_targets_neuron,
            "neuron_share": neuron_share,
            "top_5_targets": top5_targets,
            "regulating_tfs": regulating_tfs,
            "rank_position": rank_pos.get(v6_id, 0),
            "integrated_score": round(float(scores.get(v6_id, 0)), 4),
            "proof_status": proof.get(v6_id, ""),
        })

    out_df = pd.DataFrame(records)
    # Rank regulators by neuron-fate activity first (this is a neural-TF
    # project), then by unique-target out-degree.
    out_df = out_df.sort_values(
        ["n_targets_neuron", "n_targets_total"], ascending=False
    ).reset_index(drop=True)

    # ---------- Assertions ---------------------------------------------------
    assert len(out_df) > 0, "ANANSE scan produced zero rows"
    req_cols = ["v6_id", "gene_name", "h1smcg_id", "is_ananse_tf",
                "is_ananse_target", "fates_as_tf", "n_targets_total",
                "n_targets_neuron", "neuron_share", "top_5_targets",
                "regulating_tfs", "rank_position", "integrated_score",
                "proof_status"]
    for col in req_cols:
        assert col in out_df.columns, f"Missing required column: {col}"

    # ---------- Write outputs -------------------------------------------------
    RESULTS.mkdir(parents=True, exist_ok=True)

    full_path = RESULTS / "ananse_network_full.csv"
    out_df.to_csv(full_path, index=False)
    print(f"\n  Full network scan: {full_path} ({len(out_df)} candidates)")

    parquet_path = RESULTS / "ananse_network_full.parquet"
    try:
        out_df.to_parquet(parquet_path, index=False)
        print(f"  Parquet backup: {parquet_path}")
    except Exception as e:
        print(f"  [warn] Parquet write failed: {e}")

    top_regulators = (
        out_df[out_df["is_ananse_tf"]]
        .sort_values(
            ["n_targets_neuron", "n_targets_total"], ascending=False
        )
        .head(20)
        [["v6_id", "gene_name", "h1smcg_id", "fates_as_tf",
          "n_targets_total", "n_targets_neuron", "neuron_share",
          "top_5_targets", "integrated_score", "proof_status"]]
        .reset_index(drop=True)
    )
    top_reg_path = RESULTS / "ananse_top_regulators.csv"
    top_regulators.to_csv(top_reg_path, index=False)
    print(f"  Top regulators (n={len(top_regulators)}): {top_reg_path}")

    # Summary stats
    n_tf = out_df["is_ananse_tf"].sum()
    n_tgt = out_df["is_ananse_target"].sum()
    n_both = (out_df["is_ananse_tf"] & out_df["is_ananse_target"]).sum()
    n_neither = (~out_df["is_ananse_tf"] & ~out_df["is_ananse_target"]).sum()
    print(f"\n  Summary:")
    print(f"    ANANSE TF nodes         : {n_tf}/{len(out_df)}")
    print(f"    ANANSE target genes      : {n_tgt}/{len(out_df)}")
    print(f"    Both TF and target       : {n_both}/{len(out_df)}")
    print(f"    Neither (not in network) : {n_neither}/{len(out_df)}")
    print(f"    h1SMcG mapping rate      : {n_mapped}/{len(v6_ids)}")

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

