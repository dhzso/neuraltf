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


def _auc_se_hanley(auc: float, n1: int, n2: int) -> float:
    """Hanley-McNeil standard error of the AUC.

    Q1 = AUC/(2-AUC), Q2 = 2*AUC^2/(1+AUC) (Hanley & McNeil 1982)."""
    if not (0.0 < auc < 1.0):
        return float("nan")
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc * auc / (1.0 + auc)
    return float(np.sqrt(
        (auc * (1 - auc) + (n1 - 1) * (q1 - auc * auc)
         + (n2 - 1) * (q2 - auc * auc)) / (n1 * n2)
    ))


def _honest_scores(rank: pd.DataFrame) -> pd.DataFrame:
    """Add the honest score column (label-encoding streams excluded).

    2026-09-13: mirrors precision_recall.recompute_excluding_circular —
    the ortholog headline previously used integrated_score, which
    contains rnai and perez_lineage, i.e. (near-)copies of the gold
    standard's own neural-family classification.
    """
    STREAMS = ["expression", "specificity", "reproducibility", "rnai",
               "correlation", "neural_enriched", "neural_specificity",
               "perez_lineage", "perez_influence", "fincher_brain",
               "cui_temporal"]
    W = {"expression": 0.1, "specificity": 0.1, "reproducibility": 0.1,
         "rnai": 0.05, "correlation": 0.05, "neural_enriched": 0.1,
         "neural_specificity": 0.1, "perez_lineage": 0.1,
         "perez_influence": 0.1, "fincher_brain": 0.1, "cui_temporal": 0.1}
    excl = {"rnai", "neural_enriched", "neural_specificity", "perez_lineage"}
    keep = [s for s in STREAMS if s in rank.columns and s not in excl]
    S = rank[keep].to_numpy(dtype=float)
    Wv = np.array([W[s] for s in keep])
    valid = ~np.isnan(S)
    num = np.nan_to_num(S, nan=0.0) @ Wv
    den = valid.astype(float) @ Wv
    den = np.where(den > 0, den, 1.0)
    out = rank.copy()
    out["_honest_score"] = num / den
    return out


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
    # honest score (2026-09-13 circularity fix; see _honest_scores)
    eval_df = _honest_scores(eval_df)
    rank = _honest_scores(rank)

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

    # 2026-09-13 circularity fix: the previous headline AUC was computed
    # on integrated_score — the FULL 11-stream score including rnai (11/47
    # positives carry the screened label itself) and perez_lineage (34/47
    # positives carry 1.0 — nearly the same neural-family classification
    # the gold standard encodes). That measured label re-encoding, not
    # cross-species validation. The honest score (rnai/neural_*/
    # perez_lineage excluded) is now the headline; the circular value is
    # retained for transparency.
    from scipy import stats as _st

    def _midrank_auc(y: np.ndarray, s: np.ndarray) -> float:
        r = _st.rankdata(s)
        n_pos = int(y.sum())
        n_neg = int(len(y) - n_pos)
        return float((r[y == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))

    y = eval_df["label"].to_numpy()
    auc_circ = _midrank_auc(y, eval_df["integrated_score"].to_numpy())
    auc_honest = _midrank_auc(y, eval_df["_honest_score"].to_numpy())

    # Mann-Whitney on the honest score (one-sided, positives greater)
    pos_scores = eval_df.loc[eval_df["label"] == 1, "_honest_score"].to_numpy()
    neg_scores = eval_df.loc[eval_df["label"] == 0, "_honest_score"].to_numpy()
    if len(pos_scores) and len(neg_scores):
        mw = _st.mannwhitneyu(pos_scores, neg_scores, alternative="greater")
        mw_u, mw_p = float(mw.statistic), float(mw.pvalue)
    else:
        mw_u, mw_p = float("nan"), float("nan")

    # Hanley-McNeil 95% CI on the honest AUC (SE per Hanley & McNeil 1982)
    n1, n2 = int(y.sum()), int(len(y) - y.sum())
    if n1 > 0 and n2 > 0:
        q1 = eval_df.loc[eval_df["label"] == 1, "_honest_score"]
        q2 = eval_df.loc[eval_df["label"] == 0, "_honest_score"]
        auc_v = auc_honest
        se = _auc_se_hanley(auc_v, n1, n2)
        ci = (max(0.0, auc_v - 1.96 * se), min(1.0, auc_v + 1.96 * se))
    else:
        se, ci = float("nan"), (float("nan"), float("nan"))

    print(f"  ROC-AUC (honest score, label streams excluded): {auc_honest:.3f} "
          f"95% CI {ci[0]:.3f}-{ci[1]:.3f}")
    print(f"  ROC-AUC (circular integrated score, transparency only): {auc_circ:.3f}")
    print(f"  Mann-Whitney U (honest, one-sided): {mw_u:.1f}  p={mw_p:.3e}")

    # top-decile enrichment of positives (honest score), WITH the
    # negative rate shown (previously suppressed — the printed
    # positives-only rate read like enrichment evidence while negatives
    # were also enriched).
    q = rank["_honest_score"].quantile(0.9)
    top = set(rank.loc[rank["_honest_score"] >= q, "gene_id"])
    pos_in = pos & top
    neg_in = neg & top
    pos_cov = len(pos & set(rank["gene_id"]))
    neg_cov = len(neg & set(rank["gene_id"]))
    pos_rate = len(pos_in) / max(pos_cov, 1)
    neg_rate = len(neg_in) / max(neg_cov, 1)
    print(f"  top-decile rate (honest): positives {len(pos_in)}/{pos_cov} "
          f"({pos_rate:.2f}) vs negatives {len(neg_in)}/{neg_cov} "
          f"({neg_rate:.2f})")

    result = {
        "coverage_n_positive": int(len(pos)),
        "coverage_n_negative": int(len(neg)),
        "covered_in_rank": int(len(covered)),
        "roc_auc": None if np.isnan(auc_honest) else float(auc_honest),
        "roc_auc_ci95": [None if np.isnan(ci[0]) else float(ci[0]),
                         None if np.isnan(ci[1]) else float(ci[1])],
        "roc_auc_circular_transparency": None if np.isnan(auc_circ) else float(auc_circ),
        "mann_whitney_u": None if np.isnan(mw_u) else float(mw_u),
        "mann_whitney_p_one_sided": None if np.isnan(mw_p) else float(mw_p),
        "top_decile_positive_rate": pos_rate,
        "top_decile_negative_rate": neg_rate,
        "note": (
            "Provisional, coverage-limited. Built from PlanMine Homo-sapiens BLAST "
            "descriptions only (~456 annotated genes); a full ortholog table "
            "(Ensembl Compara / DIOPT) is required for a definitive benchmark. "
            "2026-09-13: the headline AUC is computed on the honest score "
            "(rnai/neural_enriched/neural_specificity/perez_lineage excluded); "
            "the previous integrated-score AUC measured label re-encoding "
            "(34/47 positives carried perez_lineage=1.0, 11/47 carried rnai=1)."
        ),
    }
    with open(RESULTS_DIR / "ortholog_benchmark.json", "w") as f:
        json.dump(result, f, indent=2)
    print(f"  saved: {RESULTS_DIR / 'ortholog_gold_standard.csv'}")
    return 0 if not np.isnan(auc_honest) else 1


if __name__ == "__main__":
    raise SystemExit(main())