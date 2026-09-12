#!/usr/bin/env python
"""Independent cross-species ortholog benchmark for the NeuralTF ranking.

PROVISIONAL / LIMITED-COVERAGE (documented honestly):
    The ONLY cross-species ortholog signal available in-repo is the PlanMine
    BLAST annotation parquet (``datasets/processed/planmine_annotations.parquet``),
    whose ``blast`` rows cover ~456 planarian genes and use descriptive protein
    names (e.g. ``empty spiracles homeobox 2``), NOT clean gene symbols. A full
    human<->planarian ortholog table (Ensembl Compara / DIOPT / OMA) is a
    data-acquisition task for the next milestone; this script builds the best
    defensible gold standard from what is actually on disk.

Classification:
    Each Homo-sapiens BLAST description is classified by an explicit, curated
    list of name fragments into one of:

      positive  - conserved neural-fate / neural-lineage transcription factor
                  (proneural bHLH, homeodomain patterning, glial, neural-crest)
      negative  - transcription factor of a NON-neural lineage or a general
                  TF factor (endoderm/mesoderm TFs, nuclear receptors,
                  muscle/heart/blood TFs, generic zinc fingers)
      ignored   - solute carriers, enzymes, structural proteins, or ambiguous
                  factors (not transcription-factor fate evidence either way)

    Fragments are curated from the actual HS descriptions and map to established
    gene nomenclature (ATOH/NEUROD/OLIG2/PAX6/DLX/EMX/LHX/NKX2/SIX/... positive;
    HNF/MEF2/RUNX/GATA/THR/RXR/FOXA/... negative). The full classification is
    deterministic and reproducible; a gene matching BOTH a positive and a
    negative fragment is dropped (ambiguous).

Outputs:
    results/ortholog_gold_standard.csv   (v6_id, label, description)
    A 2-class ROC-AUC (positive vs negative) over the coverage-limited universe,
    plus top-decile enrichment of positives, printed and saved to
    results/ortholog_benchmark.json.

Usage:
    python scripts/stats/ortholog_benchmark.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
PLANMINE = REPO / "datasets" / "processed" / "planmine_annotations.parquet"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# --- curated description-fragment classification ---------------------------------
POSITIVE_FRAGMENTS = [
    "atonal",                          # ATOH proneural bHLH
    "neuronal differentiation",        # NEUROD
    "oligodendrocyte lineage",         # OLIG2 (glial/neural)
    "brain specific homeobox",         # BSX
    "castor zinc finger",              # CASZ1
    "cut-like homeobox",               # CUX1
    "distal-less homeobox",            # DLX
    "dorsal root ganglia homeobox",    # DRGX
    "early b-cell factor",             # EBF
    "early growth response",           # EGR (activity-dependent neural)
    "empty spiracles homeobox",        # EMX2
    "even-skipped homeobox",           # EVX1
    "glial cells missing",             # GCM (glial)
    "hes family",                      # HES
    "hes related",                     # HEY
    "iroquois homeobox",               # IRX
    "lim homeobox",                    # LHX / LMX
    "lymphoid enhancer binding",       # LEF1 (neural crest / Wnt)
    "msh homeobox",                    # MSX1 (neural crest)
    "nescient helix-loop-helix",       # NHLH1
    "nk2 homeobox",                    # NKX2
    "one cut homeobox",                # ONECUT
    "paired box 6",                    # PAX6
    "paired-like homeodomain 3",       # PITX3 (dopaminergic)
    "pancreas specific transcription factor",  # PTF1A (GABAergic)
    "prospero homeobox",               # PROX1
    "short stature homeobox",          # SHOX2
    "single-minded family",            # SIM
    "spalt-like",                      # SALL
    "six homeobox",                    # SIX
    "barh like homeobox",              # BARHL
    "pbx/knotted 1 homeobox",          # PKNOX (neural cofactor)
]

NEGATIVE_FRAGMENTS = [
    "hepatocyte nuclear factor",       # HNF (endoderm/liver)
    "myocyte enhancer factor",         # MEF2 (muscle)
    "runt-related",                    # RUNX (hematopoietic)
    "scleraxis",                       # SCX (tendon/mesoderm)
    "snail family",                    # SNAI (mesenchyme/EMT)
    "thyroid hormone receptor",        # THRA/THRB
    "retinoid x receptor",             # RXR
    "interferon regulatory factor",    # IRF (immune)
    "gata binding protein",            # GATA (heart/blood)
    "smad family",                     # SMAD (signaling)
    "paired box 5",                    # PAX5 (B-cell)
    "paired box 8",                    # PAX8 (thyroid/kidney)
    "nuclear receptor subfamily",      # NR2… (nuclear receptors)
    "homeobox b",                      # HOXB (axial, not neural-fate)
    "homeobox c",                      # HOXC (axial)
    "forkhead box a",                  # FOXA (endoderm)
    "forkhead box c",                  # FOXC (mesoderm)
    "forkhead box f",                  # FOXF (mesoderm)
    "forkhead box j",                  # FOXJ (ciliary)
    "forkhead box k",                  # FOXK (general)
    "forkhead box l",                  # FOXL (mesenchyme)
    "forkhead box o",                  # FOXO (general)
    "forkhead box d",                  # FOXD (mesoderm/neural-crest)
    "histone h4 transcription factor", # HINFP (general)
    "regulatory factor x",             # RFX (general)
    "nuclear transcription factor y",  # NFY (general)
    "ovo-like",                        # OVOL (epidermis)
    "transcription factor 4",          # TCF4 (broad)
    "transcription factor 7-like",     # TCF7L2 (broad)
    "transcription factor dp-1",       # TFDP1 (general)
    "doublesex and mab-3",             # DMRT (gonadal)
    "interferon regulatory",           # IRF
]


def classify(desc: str) -> str | None:
    d = (desc or "").lower()
    if not d:
        return None
    pos = any(p in d for p in POSITIVE_FRAGMENTS)
    neg = any(n in d for n in NEGATIVE_FRAGMENTS)
    if pos and neg:
        return None  # ambiguous
    if pos:
        return "positive"
    if neg:
        return "negative"
    return None  # ignored (not TF-fate evidence)


def build_gold_standard() -> pd.DataFrame:
    if not PLANMINE.exists():
        raise FileNotFoundError(f"{PLANMINE} missing; run query_planmine.py first")
    ann = pd.read_parquet(PLANMINE)
    blast = ann[ann["kind"] == "blast"].copy()
    hs = blast[blast["value"].astype(str).str.lower().str.contains("homo sapiens", na=False)]

    rows = []
    for gid, sub in hs.groupby("gene_id_v6"):
        descs = [str(d) for d in sub["description"].dropna().astype(str) if str(d).lower() != "nan"]
        for d in descs:
            lbl = classify(d)
            if lbl is not None:
                rows.append({"gene_id": gid, "label": lbl, "description": d})
    df = pd.DataFrame(rows).drop_duplicates(subset=["gene_id", "label"])
    # A gene could be positive via one description and negative via another:
    # collapse to one row per gene, dropping conflicting labels.
    dedup = []
    for gid, sub in df.groupby("gene_id"):
        labels = sorted(set(sub["label"]))
        if len(labels) == 1:
            dedup.append({"gene_id": gid, "label": labels[0],
                          "description": sub.iloc[0]["description"]})
    return pd.DataFrame(dedup)


def roc_auc(y_true: np.ndarray, scores: np.ndarray) -> float:
    order = np.argsort(-scores, kind="stable")
    y = y_true[order]
    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    tpr, fpr = [0.0], [0.0]
    tp = fp = 0
    for lbl in y:
        if lbl:
            tp += 1
        else:
            fp += 1
        tpr.append(tp / n_pos)
        fpr.append(fp / n_neg)
    return float(np.trapezoid(tpr, fpr))


def main() -> int:
    rank = pd.read_csv(RUN_DIR / "rank.csv").drop_duplicates(subset="gene_id", keep="first")
    gs = build_gold_standard()
    gs.to_csv(RESULTS_DIR / "ortholog_gold_standard.csv", index=False)

    pos = set(gs[gs["label"] == "positive"]["gene_id"])
    neg = set(gs[gs["label"] == "negative"]["gene_id"])

    # Restrict to covered genes present in the rank table
    covered = sorted((pos | neg) & set(rank["gene_id"]))
    eval_df = rank[rank["gene_id"].isin(covered)].copy()
    eval_df["label"] = eval_df["gene_id"].isin(pos).astype(int)

    tested_pos = pos & set(rank[rank["proof_status"] == "tested"]["gene_id"])
    conf_pos = pos & set(rank.loc[
        rank.get("phenotype_confirmed", pd.Series(False, index=rank.index))
        .astype(bool), "gene_id"])
    print("Cross-species ortholog benchmark (Homo sapiens BLAST descriptions)")
    print(f"  positive (neural-fate orthologs):  {len(pos)} labeled, {len(pos & set(rank['gene_id']))} in rank")
    print(f"  negative (non-neural orthologs):   {len(neg)} labeled, {len(neg & set(rank['gene_id']))} in rank")
    print(f"  covered (in-rank, classified):     {len(covered)}")
    print(f"  tested-gene overlap among positives (shared biology, not leakage): "
          f"{len(tested_pos)}")
    print(f"  phenotype-confirmed overlap among positives: {len(conf_pos)}")

    auc = roc_auc(eval_df["label"].to_numpy(), eval_df["integrated_score"].to_numpy())
    print(f"  ROC-AUC (positive vs negative, integrated score): {auc:.3f}")

    # top-decile enrichment of positives
    q = rank["integrated_score"].quantile(0.9)
    top = rank[rank["integrated_score"] >= q]["gene_id"]
    n_pos_top = len(pos & set(top))
    prevalence = len(pos & set(rank["gene_id"])) / max(len(rank), 1)
    print(f"  positives in top decile: {n_pos_top}/{len(pos & set(rank['gene_id']))} "
          f"(prevalence {prevalence:.4f})")

    result = {
        "coverage_n_positive": int(len(pos)),
        "coverage_n_negative": int(len(neg)),
        "covered_in_rank": int(len(covered)),
        "roc_auc": None if np.isnan(auc) else float(auc),
        "note": (
            "Provisional, coverage-limited. Built from PlanMine Homo-sapiens BLAST "
            "descriptions only (~456 annotated genes); a full ortholog table "
            "(Ensembl Compara / DIOPT) is required for a definitive benchmark."
        ),
    }
    with open(RESULTS_DIR / "ortholog_benchmark.json", "w") as f:
        json.dump(result, f, indent=2)
    print(f"  saved: {RESULTS_DIR / 'ortholog_gold_standard.csv'}")
    return 0 if not np.isnan(auc) else 1


if __name__ == "__main__":
    raise SystemExit(main())