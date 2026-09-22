#!/usr/bin/env python
"""Annotate all NeuralTF outputs with the corrected King-2024 ground truth.

2026-09-11 ground-truth correction (review follow-up): King mmc5 lists ALL
TFs inhibited by RNAi ("All Transcription Factors Inhibited"); the
distributed xlsx copy lost the red/green phenotype font encoding, so the
historical "tested" label meant SCREENED, not phenotype-validated. This
script re-derives the FISH-confirmed subset from the paper's figures
(bioforge.evidence.groundtruth.PHENOTYPE_CONFIRMED_V6) and stamps:

    phenotype_confirmed (bool)  — FISH-confirmed loss-of-cell-type phenotype
    ground_truth_status (str)    — phenotype_confirmed | screened_no_phenotype |
                                   not_screened

onto every result table WITHOUT re-running the expensive pipeline (stream
scores are unchanged; only labels/annotations are corrected). Idempotent:
re-running overwrites the columns with the same values.

Outputs touched (all under projects/NeuralTF):
    runs/pipeline_run/rank.csv               (+phenotype_confirmed)
    runs/pipeline_run/rank_neural.csv       (+phenotype_confirmed, ground_truth_status)
    runs/pipeline_run/pipeline_results.json (+phenotype_confirmed per candidate)
    runs/pipeline_run/evidence_cards.md      (proof-status lines re-rendered)
    results/fixed_full_rank.csv              (+phenotype_confirmed)
    results/dirichlet_centered_full_rank.csv (+phenotype_confirmed)
    results/dirichlet_uniform_full_rank.csv  (+phenotype_confirmed)
    results/top10_neural_tfs_prioritized.csv (+phenotype_confirmed, notes fix)
    results/dirichlet_centered_top10.csv     (+phenotype_confirmed, notes fix)
    results/dirichlet_uniform_top10.csv      (+phenotype_confirmed, notes fix)
    results/dirichlet_centered_overall_top10.csv / uniform_overall_top10.csv
    results/candidate_summary_report.md      (wording + flags)
    results/supplementary_table_S1..S4       (+phenotype_confirmed)
    results/tf_ranked_*.csv                  (+phenotype_confirmed where gene ids present)

Usage:
    python projects/NeuralTF/scripts/annotate_phenotype_groundtruth.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from bioforge.evidence.groundtruth import (  # noqa: E402
    MMC5_GROUND_TRUTH_NOTE,
    PHENOTYPE_CONFIRMED_V6,
    phenotype_status,
)

RUN = ROOT / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RES = ROOT / "projects" / "NeuralTF" / "results"


def _stamp(df: pd.DataFrame, gene_col: str, add_status: bool = True) -> pd.DataFrame:
    """Add phenotype_confirmed (and optionally ground_truth_status) columns."""
    out = df.copy()
    gid = out[gene_col].astype(str).str.strip()
    confirmed = gid.isin(PHENOTYPE_CONFIRMED_V6)
    out["phenotype_confirmed"] = confirmed
    if add_status:
        rnai = pd.to_numeric(out["rnai"], errors="coerce") \
            if "rnai" in out.columns else pd.Series(0.0, index=out.index)
        out["ground_truth_status"] = [
            phenotype_status(r, g) for r, g in zip(rnai.fillna(0.0), gid)
        ]
    return out


def _fix_shortlist_csv(path: Path) -> bool:
    """Shortlist CSVs (fixed/centered/uniform top10): add flags and fix the
    'Known validated neural TF' note to distinguish confirmed vs screened."""
    if not path.exists():
        return False
    df = pd.read_csv(path)
    gene_col = "gene_id_v6" if "gene_id_v6" in df.columns else (
        "gene_id" if "gene_id" in df.columns else None)
    if gene_col is None:
        print(f"  [skip] {path.name}: no gene id column")
        return False
    df = _stamp(df, gene_col, add_status=False)
    note_col = "rnai_screen_or_marker_notes"
    if note_col in df.columns:
        def _rewrite(r):
            s = str(r.get(note_col, ""))
            if "Known validated neural TF" in s:
                base = ("FISH-confirmed phenotype TF (King 2024)"
                        if bool(r.get("phenotype_confirmed"))
                        else "RNAi-screened in King 2024 (no published phenotype)")
                s = s.replace("Known validated neural TF (positive control)", base)
            return s
        df[note_col] = df.apply(_rewrite, axis=1)
    df.to_csv(path, index=False)
    print(f"  [ok] {path.name} ({len(df)} rows, "
          f"{int(df['phenotype_confirmed'].sum())} phenotype-confirmed)")
    return True


def main() -> int:
    print("=== Annotating outputs with corrected King-2024 ground truth ===")
    print(f"  Phenotype-confirmed set: {len(PHENOTYPE_CONFIRMED_V6)} genes")

    # ---- rank.csv / rank_neural.csv -----------------------------------
    for name, add_status in (("rank.csv", False), ("rank_neural.csv", True)):
        p = RUN / name
        if not p.exists():
            print(f"  [skip] {name} missing")
            continue
        df = pd.read_csv(p)
        df = _stamp(df, "gene_id", add_status=add_status)
        df.to_csv(p, index=False)
        n_conf = int(df["phenotype_confirmed"].sum())
        if "ground_truth_status" in df.columns:
            print(f"  [ok] {name}: {n_conf} phenotype-confirmed | "
                  f"{(df['ground_truth_status'] == 'screened_no_phenotype').sum()} "
                  f"screened-no-phenotype | not_screened "
                  f"{(df['ground_truth_status'] == 'not_screened').sum()}")
        else:
            print(f"  [ok] {name}: {n_conf} phenotype-confirmed")

    # ---- pipeline_results.json ----------------------------------------
    p = RUN / "pipeline_results.json"
    if p.exists():
        payload = json.loads(p.read_text())
        n_fixed = 0
        for cand in payload.get("top_candidates", []):
            cand["phenotype_confirmed"] = cand["gene_id"] in PHENOTYPE_CONFIRMED_V6
            n_fixed += int(cand["phenotype_confirmed"])
        payload["ground_truth_note"] = MMC5_GROUND_TRUTH_NOTE
        p.write_text(json.dumps(payload, indent=2))
        print(f"  [ok] pipeline_results.json ({n_fixed} confirmed in top-50)")
    else:
        print("  [skip] pipeline_results.json missing")

    # ---- evidence_cards.md: re-render proof-status lines ---------------
    # Line-based (the earlier per-gene regex with .*? + re.S over the 9 MB
    # file backtracked catastrophically; a card's Proof-status line is
    # always within a few lines of its Gene-ID line, so a bounded
    # single-pass scan is both fast and exact).
    p = RUN / "evidence_cards.md"
    if p.exists():
        text = p.read_text(encoding="utf-8")
        if "GROUND-TRUTH NOTE (2026-09-11)" not in text:
            note = ("<!-- GROUND-TRUTH NOTE (2026-09-11): 'tested' = RNAi-screened "
                    "in King 2024 mmc5 ('All Transcription Factors Inhibited'), "
                    "NOT phenotype-validated. FISH-confirmed phenotypes are marked "
                    "'(phenotype-confirmed †)' on the card. -->\n\n")
            text = note + text
        lines = text.split("\n")
        n_stamped = 0
        cur_gene: str | None = None
        for i, ln in enumerate(lines):
            m_gid = re.match(r"- \*\*Gene ID:\*\* `([^`]+)`", ln)
            if m_gid:
                cur_gene = m_gid.group(1)
                continue
            if (ln == "- **Proof status:** tested" and cur_gene
                    and cur_gene in PHENOTYPE_CONFIRMED_V6):
                lines[i] = "- **Proof status:** tested (phenotype-confirmed †)"
                n_stamped += 1
        p.write_text("\n".join(lines), encoding="utf-8")
        print(f"  [ok] evidence_cards.md re-stamped ({n_stamped} cards marked)")
    else:
        print("  [skip] evidence_cards.md missing")

    # ---- full-rank CSVs -------------------------------------------------
    for name in ("fixed_full_rank.csv", "dirichlet_centered_full_rank.csv",
                 "dirichlet_uniform_full_rank.csv"):
        p = RES / name
        if p.exists():
            df = pd.read_csv(p)
            df = _stamp(df, "gene_id", add_status=False)
            df.to_csv(p, index=False)
            print(f"  [ok] {name} "
                  f"({int(df['phenotype_confirmed'].sum())} phenotype-confirmed)")
        else:
            print(f"  [skip] {name} missing")

    # ---- shortlists ------------------------------------------------------
    for name in ("top10_neural_tfs_prioritized.csv",
                 "dirichlet_centered_top10.csv", "dirichlet_uniform_top10.csv",
                 "dirichlet_centered_overall_top10.csv",
                 "dirichlet_uniform_overall_top10.csv"):
        _fix_shortlist_csv(RES / name)

    # ---- supplementary tables -------------------------------------------
    for name in ("supplementary_table_S1_method_comparison.csv",
                 "supplementary_table_S2_fixed_all_candidates.csv",
                 "supplementary_table_S3_centered_all_candidates.csv",
                 "supplementary_table_S4_uniform_all_candidates.csv",
                 "supplementary_table_S5_tf_neural.csv",
                 "supplementary_table_S6_tf_all.csv",
                 "supplementary_table_S7_tf_catalog.csv"):
        p = RES / name
        if p.exists():
            df = pd.read_csv(p)
            gene_col = "gene_id_v6" if "gene_id_v6" in df.columns else "gene_id"
            if gene_col not in df.columns:
                gene_col = df.columns[0]  # S5-S7 use v6 ids in col 0
            df = _stamp(df, gene_col, add_status=False)
            df.to_csv(p, index=False)
            print(f"  [ok] {name}")
        else:
            print(f"  [skip] {name} missing")

    # ---- TF ranked tables --------------------------------------------------
    for name in ("tf_ranked_neural_top19.csv", "tf_ranked_all_top43.csv",
                 "tf_ranked_catalog_top74.csv"):
        p = RES / name
        if p.exists():
            df = pd.read_csv(p)
            cand = [c for c in ("gene_id_v6", "gene_id", "v6_id") if c in df.columns]
            if cand:
                df = _stamp(df, cand[0], add_status=False)
                df.to_csv(p, index=False)
                print(f"  [ok] {name}")
            else:
                print(f"  [skip] {name}: no gene column")
        else:
            print(f"  [skip] {name} missing")

    print("\nDone. NOTE: " + MMC5_GROUND_TRUTH_NOTE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
