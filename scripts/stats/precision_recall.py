#!/usr/bin/env python
"""Precision-recall analysis with multiple ground-truth labels.

WS3 fix for circularity: the `rnai` and `neural_enriched` streams ARE
part of the integrated score, and the primary ground truth (proof_status ==
tested) is derived from the same King mmc5 RNAi table that feeds the rnai
stream — so a naive evaluation is circular and inflates ROC-AUC
(~0.91 vs ~0.69 honest). This script therefore reports:

  1. circular   — score includes rnai + neural_enriched streams,
                  label = "tested" (RNAi-screened in King 2024 mmc5)
  2. honest     — score recomputed EXCLUDING rnai, neural_enriched,
                  neural_specificity, perez_lineage (streams directly
                  encoding the labels)
  3. strict     — honest additionally excludes reproducibility
                  (atlas membership embeds King-neural-G0 hits)

Each of the above is evaluated against BOTH labels (2026-09-11
ground-truth correction):
  - "screened"            proof_status == tested (King 2024 mmc5 RNAi
                          screening list; phenotype NOT implied — mmc5 is
                          titled "All Transcription Factors Inhibited"
                          and the distributed copy lost the red/green
                          phenotype font encoding)
  - "phenotype_confirmed" FISH-confirmed loss-of-cell-type phenotypes
                          (paper Fig 3J/4E, S4, S7, S8; the strictest
                          defensible "validated" label), via
                          bioforge.evidence.groundtruth

Also reports a per-stream AUC table (both labels) and a
presence-conditioned AUC diagnostic for fincher_brain — its
missing-vs-present pattern carries signal that an unconditional
fillna(0) AUC hides.

The phenotype_confirmed numbers are the publishable "recovery of
validated neural TFs" estimates; the screened numbers describe recovery
of the King screening list.

Outputs:
   results/precision_recall.json (curves + per-stream AUCs, both labels)

2026-09-19 tie fix: ROC/PR curves and AP are computed with GROUPED
THRESHOLDS — all candidates sharing a score form one tie-averaged point
(previously the curve stepped through tied rows in arbitrary file order,
making the trapezoid AUCs row-order dependent given the documented
massive tie groups). Exported curve arrays are now one point per
distinct score threshold. ROC-AUC is now exactly the tie-correct
midrank (Mann-Whitney) AUC.

2026-09-26: precision@k is tie-aware as well — when the k-th row falls
inside a tie group, the reported value is the expectation over a
uniformly random ordering of that group (hypergeometric mean of the
positives inside the cut), matching the grouped-threshold convention
instead of a file-order slice.

Usage:
    python scripts/stats/precision_recall.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from bioforge.evidence.groundtruth import (  # noqa: E402
    MMC5_GROUND_TRUTH_NOTE,
    PHENOTYPE_CONFIRMED_V6,
)
from bioforge.evidence.scoring import DEFAULT_WEIGHTS, STREAM_ORDER  # noqa: E402

RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

STREAMS = [getattr(s, "value", s) for s in STREAM_ORDER]
W_DEFAULT = {getattr(k, "value", k): float(v)
             for k, v in DEFAULT_WEIGHTS.items()}
# streams that directly encode the ground-truth labels. 2026-09-06 audit:
# perez_lineage is a hand-curated neural-TF-family membership list (AUC
# 0.9903 alone on the RNAi-validated label — the same species of leakage
# as rnai) and MUST be excluded from the honest score.
CIRCULAR_STREAMS = {"rnai", "neural_enriched", "neural_specificity",
                    "perez_lineage"}
# Stricter honest variant: additionally drops reproducibility (atlas
# membership embeds King-neural-G0 hits; AUC 0.9078 alone). Reported as
# the lower bound on true discrimination.
CIRCULAR_STREAMS_STRICT = CIRCULAR_STREAMS | {"reproducibility"}
# 2026-09-13 king_free variant: excludes EVERY stream that draws on King
# 2024 data (mmc5/6/7). The audit found the honest/strict AUCs still
# retain King-study evidence — binary King-atlas membership ALONE scores
# AUC ~0.837 (screened) / ~0.808 (phenotype_confirmed), nearly the full
# honest 0.846/0.875, because expression/specificity carry King-mmc7
# floors and correlation is King mmc6. king_free keeps only streams
# derived from OTHER studies:
#   expression/specificity still contain King floors (they are max-fused
#   per stream, so a per-stream exclusion is impossible without a
#   pipeline re-run; noted as a residual caveat) — the reported king_free
#   number is therefore still an UPPER bound, but with the pure King
#   streams (rnai, correlation, neural_enriched, neural_specificity) and
#   the King-membership reproducibility term removed, it isolates how
#   much discrimination the non-King studies (Fincher, Plass, Cui,
#   Perez-influence, Fincher-brain) provide on their own.
KING_FREE_EXCLUDE = CIRCULAR_STREAMS_STRICT | {"correlation"}


def load_candidates():
    p = RUN_DIR / "rank.csv"
    if not p.exists():
        raise FileNotFoundError(f"No candidate score file at {p}")
    return pd.read_csv(p).drop_duplicates(subset="gene_id", keep="first")


def recompute_excluding_circular(df: pd.DataFrame,
                                 exclude: set[str] | None = None) -> np.ndarray:
    """Renormalized weighted score over the non-circular streams only.

    ``exclude`` defaults to CIRCULAR_STREAMS; pass
    CIRCULAR_STREAMS_STRICT for the reproducibility-free variant.
    """
    excl = CIRCULAR_STREAMS if exclude is None else exclude
    keep = [s for s in STREAMS if s in df.columns and s not in excl]
    S = df[keep].to_numpy(dtype=float)
    W = np.array([W_DEFAULT[s] for s in keep])
    valid = ~np.isnan(S)
    S_filled = np.nan_to_num(S, nan=0.0)
    num = S_filled @ W
    den = valid.astype(float) @ W
    den = np.where(den > 0, den, 1.0)
    return num / den


def _group_ends(y_scores_sorted: np.ndarray) -> np.ndarray:
    """End-indices (exclusive) of each distinct-score threshold group.

    Scores must already be sorted (descending). Tied scores form ONE
    threshold group so the curves are tie-averaged instead of following
    an arbitrary row-order staircase inside a tie group (the 2026-09-13
    audit documented 541 genes tied at 0.5667 and 1,817 genes in tie
    groups >= 5 — the staircase made the trapezoid AUCs row-order
    dependent).
    """
    n = len(y_scores_sorted)
    if n == 0:
        return np.array([], dtype=int)
    change = np.nonzero(np.diff(y_scores_sorted))[0] + 1
    return np.concatenate([change, [n]])


def compute_pr_curve(y_true, y_scores):
    """Tie-averaged precision-recall curve (one point per threshold group)."""
    y_true = np.asarray(y_true)
    y_scores = np.asarray(y_scores)
    order = np.argsort(-y_scores, kind="stable")
    y_sorted = y_true[order]
    s_sorted = y_scores[order]
    n_pos = float(y_true.sum())
    if n_pos == 0:
        return np.array([0.0]), np.array([1.0])
    ends = _group_ends(s_sorted)
    tp_cum = np.asarray(np.cumsum(y_sorted))
    recalls = tp_cum[ends - 1] / n_pos
    precisions = tp_cum[ends - 1] / ends
    return recalls, precisions


def _trapezoid(y, x):
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    if hasattr(np, "trapz"):
        return float(np.trapz(y, x))
    y = np.asanyarray(y); x = np.asanyarray(x)
    return float(np.sum((x[1:] - x[:-1]) * (y[1:] + y[:-1]) / 2.0))


def compute_roc_curve(y_true, y_scores):
    """Tie-averaged ROC curve (one point per threshold group, from (0,0))."""
    y_true = np.asarray(y_true)
    y_scores = np.asarray(y_scores)
    order = np.argsort(-y_scores, kind="stable")
    y_sorted = y_true[order]
    s_sorted = y_scores[order]
    n_pos = int(y_true.sum())
    n_neg = int(len(y_true) - n_pos)
    ends = _group_ends(s_sorted)
    tp_cum = np.asarray(np.cumsum(y_sorted))
    fp_cum = np.asarray(np.cumsum(1 - y_sorted))
    tpr_pts = np.concatenate([[0.0], tp_cum[ends - 1] / (n_pos if n_pos else 1)])
    fpr_pts = np.concatenate([[0.0], fp_cum[ends - 1] / (n_neg if n_neg else 1)])
    return fpr_pts, tpr_pts


def _average_precision(y_true, y_scores) -> float:
    """sklearn-style AP: step-wise sum over threshold groups.

    AP = sum over groups g containing positives of
    (recall_g - recall_{g-1}) * precision_g, with precision/recall taken
    at group END (tie-averaged), not at arbitrary row order within ties.
    """
    y_true = np.asarray(y_true)
    y_scores = np.asarray(y_scores)
    order = np.argsort(-y_scores, kind="stable")
    y_sorted = y_true[order]
    s_sorted = y_scores[order]
    n_pos = float(y_true.sum())
    if n_pos == 0:
        return 0.0
    ends = _group_ends(s_sorted)
    tp_cum = np.asarray(np.cumsum(y_sorted))
    rec = tp_cum[ends - 1] / n_pos
    prec = tp_cum[ends - 1] / ends
    rec_prev = np.concatenate([[0.0], rec[:-1]])
    # only groups that add positives contribute
    d_rec = rec - rec_prev
    return float(np.sum(d_rec * prec))


def evaluate(df, score_col, label, y_true):
    n_pos = int(y_true.sum())
    y_scores = pd.to_numeric(df[score_col], errors="coerce").fillna(0).values

    recalls, precisions = compute_pr_curve(y_true, y_scores)
    fpr, tpr = compute_roc_curve(y_true, y_scores)
    ap = _average_precision(y_true, y_scores)
    pr_auc = _trapezoid(precisions, recalls) if len(recalls) > 1 else 0.0
    roc_auc = _trapezoid(tpr, fpr) if len(fpr) > 1 else 0.5

    # 2026-09-26 tie-aware precision@k: slicing the first k rows of a
    # stable sort cuts INSIDE tie groups in arbitrary file order (the
    # documented massive tie groups make precision@{5..20} row-order
    # dependent). Report the expectation over a uniformly random ordering
    # of each tie group — positives inside the cut follow the
    # hypergeometric mean — consistent with the grouped-threshold
    # curves/AP above. With no cut tie this reduces exactly to the plain
    # top-k precision.
    prec_at_k = {}
    order = np.argsort(-y_scores, kind="stable")
    y_sorted = np.asarray(y_true)[order]
    s_sorted = y_scores[order]
    ends = np.flatnonzero(np.diff(s_sorted) != 0) + 1      # tie-group ends
    ends = np.concatenate([ends, [len(s_sorted)]])
    starts = np.concatenate([[0], ends[:-1]])
    cum_tp = np.cumsum(y_sorted)
    for k in (5, 10, 15, 20):
        if k <= len(y_true):
            gi = int(np.searchsorted(ends, k - 1, side="right"))
            a, b = int(starts[gi]), int(ends[gi])
            tp_before = int(cum_tp[a - 1]) if a > 0 else 0
            tp_group = int(cum_tp[b - 1]) - tp_before
            exp_tp = tp_before + (k - a) * tp_group / (b - a)
            prec_at_k[f"precision@{k}"] = float(exp_tp / k)

    print(f"\n[{label}] n={len(y_true)}, positives={n_pos}")
    print(f"  PR-AUC: {pr_auc:.4f}  ROC-AUC: {roc_auc:.4f}  AP: {ap:.4f}")
    for k, v in prec_at_k.items():
        print(f"  {k}: {v:.4f}")

    return {
        "n_candidates": int(len(y_true)),
        "n_positives": n_pos,
        "average_precision": ap,
        "pr_auc": float(pr_auc),
        "roc_auc": float(roc_auc),
        "precision_at_k": prec_at_k,
        "score_column": score_col,
        "roc": {
            "fpr": [float(x) for x in fpr],
            "tpr": [float(x) for x in tpr],
            "auc": float(roc_auc),
        },
        "pr": {
            "recall": [float(x) for x in recalls],
            "precision": [float(x) for x in precisions],
            "auc": float(pr_auc),
            "baseline": float(n_pos / len(y_true)) if len(y_true) else 0.0,
        },
    }


def per_stream_auc(df, y_by_label):
    """AUC of each single stream against each label (leakage diagnosis).

    Also reports a presence-conditioned AUC per stream: among genes where
    the stream is PRESENT, does its value still rank positives higher?
    An unconditional fillna(0) AUC conflates presence with value when a
    stream is mostly absent for positives (e.g. fincher_brain scored
    ~0.41 unconditionally but ~0.66 present-only on the screened label).
    """
    from scipy import stats as _st

    def _auc(y, s):
        r = _st.rankdata(np.nan_to_num(s, nan=0.0))
        n_pos = int(y.sum())
        n_neg = int(len(y) - n_pos)
        if n_pos == 0 or n_neg == 0:
            return float("nan")
        return float((r[y == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))

    out = {}
    for stream in STREAMS:
        if stream not in df.columns:
            continue
        v = pd.to_numeric(df[stream], errors="coerce")
        row = {"n_present": int(v.notna().sum())}
        for lname, y in y_by_label.items():
            yv = y.astype(int).to_numpy()
            row[f"auc_{lname}"] = _auc(yv, v.to_numpy())
            m = v.notna().to_numpy()
            if m.sum() > 10 and yv[m].sum() > 0 and (yv[m] == 0).sum() > 0:
                row[f"auc_{lname}_present_only"] = _auc(yv[m], v.to_numpy()[m])
            else:
                row[f"auc_{lname}_present_only"] = float("nan")
        out[stream] = row
    return out


def main():
    print("=== Precision-Recall Analysis (circularity-controlled) ===")

    df = load_candidates()
    print(f"Loaded {len(df)} unique candidates")

    if "proof_status" not in df.columns:
        print("Error: proof_status column not found")
        return 1

    # 2026-09-11 ground-truth correction: evaluate every score variant
    # against BOTH labels.
    #   screened            = proof_status == tested (King mmc5 RNAi
    #                         screening list; phenotype NOT implied)
    #   phenotype_confirmed = FISH-confirmed loss-of-cell-type phenotype
    #                         (paper Fig 3J/4E, S4, S7, S8)
    y_screened = (df["proof_status"] == "tested").astype(int)
    y_conf = df["gene_id"].astype(str).isin(PHENOTYPE_CONFIRMED_V6).astype(int)
    labels = {
        "screened": y_screened,
        "phenotype_confirmed": y_conf,
    }

    # 1) Circular evaluation (score includes the label-bearing streams)
    circular = {
        lname: evaluate(df, "integrated_score",
                        f"circular [label={lname}] (all 11 streams)", y)
        for lname, y in labels.items()
    }

    # 2) Honest evaluation (label-encoding streams excluded — including
    #    perez_lineage, which alone carries AUC 0.99 on the label)
    df["_honest_score"] = recompute_excluding_circular(df)
    honest = {
        lname: evaluate(df, "_honest_score",
                        f"honest [label={lname}] "
                        f"(rnai/neural_*/perez_lineage excluded)", y)
        for lname, y in labels.items()
    }

    # 3) Strict honest evaluation (additionally excludes reproducibility:
    #    King-neural-G0 membership leakage via atlas_membership)
    df["_strict_score"] = recompute_excluding_circular(df, CIRCULAR_STREAMS_STRICT)
    strict = {
        lname: evaluate(df, "_strict_score",
                        f"strict [label={lname}] "
                        f"(also reproducibility excluded)", y)
        for lname, y in labels.items()
    }

    # 3b) King-free evaluation (2026-09-13): excludes ALL King-2024
    #     streams (rnai, correlation, neural_enriched, neural_specificity,
    #     + reproducibility). This is the study-of-origin control: the
    #     ground truth is King 2024, and the honest/strict scores still
    #     contained King evidence. Residual caveat: expression/specificity
    #     are max-fused with King-mmc7 floors at pipeline time, so
    #     king_free is still an upper bound on study-independent
    #     discrimination — but it isolates the non-King-study signal.
    df["_king_free_score"] = recompute_excluding_circular(df, KING_FREE_EXCLUDE)
    king_free = {
        lname: evaluate(df, "_king_free_score",
                        f"king_free [label={lname}] "
                        f"(all King-2024 streams excluded)", y)
        for lname, y in labels.items()
    }

    # 4) Per-stream leakage diagnosis (both labels, present-only variants)
    stream_aucs = per_stream_auc(df, labels)

    results = {
        "circular": circular["screened"],
        "honest": honest["screened"],
        "honest_strict": strict["screened"],
        # Keys consumers already read stay bound to the screened label for
        # continuity; the phenotype_confirmed arms are the publishable
        # "validated TF recovery" numbers.
        "phenotype_confirmed": {
            "circular": circular["phenotype_confirmed"],
            "honest": honest["phenotype_confirmed"],
            "honest_strict": strict["phenotype_confirmed"],
            "king_free": king_free["phenotype_confirmed"],
        },
        "king_free": king_free["screened"],
        "per_stream_auc": stream_aucs,
        "labels_note": (
            "Two labels are evaluated (2026-09-11): 'screened' = "
            "proof_status=='tested' (King mmc5 'All Transcription Factors "
            "Inhibited' — RNAi performed, phenotype NOT implied; the "
            "distributed mmc5 lost its red/green phenotype font encoding). "
            "'phenotype_confirmed' = FISH-confirmed loss-of-cell-type "
            "phenotypes from the paper (Fig 3J/4E, S4, S7, S8): "
            f"{len(PHENOTYPE_CONFIRMED_V6)} genes. The circular/honest/"
            "honest_strict top-level keys remain on the 'screened' label "
            "for figure compatibility; use results.phenotype_confirmed.* "
            "for any 'validated recovery' claim."
        ),
        "circularity_note": (
            "The 'circular' evaluation includes the rnai/neural_enriched/"
            "neural_specificity streams that share the King mmc5 ground "
            "truth. The 'honest' evaluation additionally excludes "
            "perez_lineage (a hand-curated neural-TF-family membership "
            "list that alone scores AUC 0.99 against the label — added in "
            "the 2026-09-06 audit after it was found leaking into every "
            "previous 'honest' contrast). The 'honest_strict' evaluation "
            "additionally excludes reproducibility (King-neural-G0 "
            "membership embedded via atlas_membership) and is the lower "
            "bound on true discrimination. Only honest/honest_strict are "
            "publishable. The 'king_free' evaluation (2026-09-13) removes "
            "every stream that draws on King 2024 (rnai, correlation, "
            "neural_enriched, neural_specificity, reproducibility) — the "
            "study-of-origin control. Residual caveat: expression/"
            "specificity are max-fused with King-mmc7 floors at pipeline "
            "time, so even king_free remains an upper bound; it isolates "
            "how much discrimination the independent studies (Fincher, "
            "Plass, Cui, Perez-influence, Fincher-brain) provide alone."
        ),
        "ground_truth_note": MMC5_GROUND_TRUTH_NOTE,
    }

    out_path = RESULTS_DIR / "precision_recall.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
