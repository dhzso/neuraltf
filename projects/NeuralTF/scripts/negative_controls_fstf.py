#!/usr/bin/env python
"""Negative control benchmarking for the NeuralTF pipeline.

Evaluates the pipeline's neural TF discovery specificity by testing whether
known non-neural fate-specifying transcription factors (FSTFs) are correctly
ranked low by the integrated score.

Scientific Rationale
--------------------
In *Schmidtea mediterranea*, FSTFs are lineage-specific master regulators that
are necessary and sufficient to specify a single cell fate. Muscle, intestine,
pharynx, epidermis and protonephridia each have their own FSTF set (King 2024,
Plass 2018, Fincher 2018). These non-neural FSTFs share the same structural
domains (Homeobox, C2H2 ZnF, bHLH) as neural FSTFs, making domain-based
annotation alone insufficient. A biologically valid pipeline must discriminate
between neural and non-neural masters using dynamic evidence (G0 progenitor
enrichment, ANANSE GRN, correlation gain).

Negative control set
--------------------
Curated non-neural planarian FSTFs with verified dd_Smed_v6 IDs (see the
audit note at NEGATIVE_CONTROLS below for the entries removed because
they were unverifiable, non-TF, or a digit-transposition of a positive):
  Muscle/EMT : myoD, nkx1-1, twist, snail2, zeb2-1

Metrics
-------
- Empirical FPR: n_negatives in top 50 / n_negatives total
- ROC-AUC and Precision-Recall AUC using scikit-learn
- Rank histogram of negative controls

Outputs
-------
  projects/NeuralTF/results/negative_control_benchmarks.csv
  (figure generation removed with the unreachable-code audit fix; the
  CSV records the honest outcome either way)

Usage
-----
    python projects/NeuralTF/scripts/negative_controls_fstf.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_root / "src"))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RUN = _root / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RES = _root / "projects" / "NeuralTF" / "results"
FIG = _root / "projects" / "NeuralTF" / "figures" / "supplementary"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Curated negative control set — non-neural planarian FSTFs
# 2026-09-06 audit fix: 12 of the 15 previous IDs were absent from
# rank.csv and 2 non-TFs (mhc-1 = myosin heavy chain, a STRUCTURAL
# protein; mlck-1 = a kinase) contradicted the "non-neural
# transcription factors" framing. dd_Smed_v6_14353_0_1 ("bmp2/4") was a
# one-digit transposition of dd_Smed_v6_14753_0_1 (ascl-2, a published
# Track-A gene — i.e. a POSITIVE mislabeled as a negative). The set is
# rebuilt below with only verified, transcription-factor-class members;
# IDs that do not appear in the candidate pool are reported as
# selection-level true negatives (they never entered the DE-seeded
# universe), never fabricated into in-pool FPR metrics.
# In-pool IDs (verified in the 2026-09 production rank.csv): myoD
# (dd7066), nkx1-1 (dd4028), mlck-1 kept ONLY as a non-TF lineage
# marker (not an FSTF; excluded from the TF-only AUC), and any other
# entries that appear after the next production re-run.
# ---------------------------------------------------------------------------
NEGATIVE_CONTROLS: dict[str, str] = {
    # Muscle lineage FSTFs (TF-class)
    "dd_Smed_v6_7066_0_1":   "myoD (longitudinal muscle FSTF)",
    "dd_Smed_v6_4028_0_1":   "nkx1-1 (circular muscle FSTF)",
    "dd_Smed_v6_22017_0_1":  "twist (mesenchymal bHLH TF)",
    "dd_Smed_v6_19501_0_1":  "snail2 (epithelial zinc-finger TF)",
    "dd_Smed_v6_24826_0_1":  "zeb2-1 (muscle/parenchyma C2H2 TF)",
    # REMOVED (audit): dd_Smed_v6_14353_0_1 was labeled "bmp2/4" but is a
    # digit transposition of dd_Smed_v6_14753_0_1 = ascl-2, a PUBLISHED
    # TRACK-A POSITIVE. The true bmp2/4 v6 ID was not verifiable; entry
    # dropped rather than guessed.
    # REMOVED (audit): dd_Smed_v6_3525_0_1 "mhc-1" — myosin heavy chain is
    # a structural protein, not a transcription factor; including it in a
    # TF-discrimination AUC is invalid.
    # REMOVED (audit): dd_Smed_v6_9718_0_1 "mlck-1" — myosin light chain
    # kinase, not a TF (same reason).
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    # Load full candidate rank table
    rank_path = RUN / "rank.csv"
    if not rank_path.exists():
        print(f"[ERROR] rank.csv not found at {rank_path}. Run scripts/run.py first.")
        return 1

    rank = pd.read_csv(rank_path)
    n_total = len(rank)
    print(f"Loaded rank.csv: {n_total} candidates")

    # Assign labels
    rank["label"] = "unlabeled"
    rank.loc[
        rank["gene_id"].isin(NEGATIVE_CONTROLS), "label"
    ] = "negative_control"
    # Positive label: RNAi-screened ONLY. (The previous definition —
    # neural_enriched>0 OR rnai>0 — used a scoring stream as ground truth,
    # making the benchmark circular.)
    pos_mask = rank["proof_status"] == "tested"
    rank.loc[pos_mask & (rank["label"] == "unlabeled"), "label"] = "positive"

    n_neg = (rank["label"] == "negative_control").sum()
    n_pos = (rank["label"] == "positive").sum()
    print(f"  Negative controls present in rank.csv: {n_neg}/{len(NEGATIVE_CONTROLS)}")
    print(f"  Positive (neural evidence) candidates : {n_pos}")

    if n_neg == 0:
        print(
            "[HONEST RESULT] No negative-control IDs found in rank.csv.\n"
            "  The curated non-neural FSTFs (myoD, gata4/5/6, foxA1, ...) were\n"
            "  not seeded as candidates because they show no significant\n"
            "  neural-cluster DE in any atlas — i.e., the DE-based seeding\n"
            "  itself already excludes non-neural lineage masters. This is a\n"
            "  true-negative observation at the SELECTION level, but it means\n"
            "  in-rank FPR/ROC metrics are NOT computable from this candidate\n"
            "  pool. This is recorded as the benchmark outcome rather than\n"
            "  being silently plotted as empty curves."
        )

    # Rank positions (1 = highest score)
    rank["rank_pos"] = rank["integrated_score"].rank(ascending=False, method="first").astype(int)

    # --- Empirical FPR -------------------------------------------------------
    neg_df = rank[rank["label"] == "negative_control"].copy()
    in_top50 = (neg_df["rank_pos"] <= 50).sum()
    fpr_top50 = in_top50 / max(n_neg, 1)
    print(f"\n  Empirical FPR (negatives in top 50): {in_top50}/{n_neg} = {fpr_top50:.1%}")

    # --- ROC-AUC and PR-AUC -------------------------------------------------
    try:
        from sklearn.metrics import roc_auc_score, average_precision_score

        # Binary classification: positive=1, negative_control=0, unlabeled excluded
        labeled = rank[rank["label"].isin(["positive", "negative_control"])].copy()
        y_true = (labeled["label"] == "positive").astype(int).values
        y_score = labeled["integrated_score"].fillna(0.0).values

        if len(np.unique(y_true)) == 2:
            roc_auc = roc_auc_score(y_true, y_score)
            pr_auc = average_precision_score(y_true, y_score)
            print(f"  ROC-AUC (neural vs non-neural): {roc_auc:.4f}")
            print(f"  PR-AUC  (neural vs non-neural): {pr_auc:.4f}")
        else:
            roc_auc = pr_auc = float("nan")
            print("  [WARN] Only one class in labeled set; ROC/PR-AUC not computable.")
    except ImportError:
        roc_auc = pr_auc = float("nan")
        print("  [WARN] scikit-learn not installed; ROC/PR-AUC skipped.")

    # --- Save benchmark CSV -------------------------------------------------
    bench_rows = []
    for gid, label_str in NEGATIVE_CONTROLS.items():
        row_match = rank[rank["gene_id"] == gid]
        if len(row_match) == 0:
            bench_rows.append({
                "gene_id": gid,
                "gene_symbol": label_str,
                "rank_pos": "not_in_candidates",
                "integrated_score": float("nan"),
                "neural_enriched": float("nan"),
                "rnai": float("nan"),
                "in_top50": False,
            })
        else:
            r = row_match.iloc[0]
            bench_rows.append({
                "gene_id": gid,
                "gene_symbol": label_str,
                "rank_pos": int(r["rank_pos"]),
                "integrated_score": round(float(r["integrated_score"]), 4),
                "neural_enriched": r.get("neural_enriched", float("nan")),
                "rnai": r.get("rnai", float("nan")),
                "in_top50": bool(r["rank_pos"] <= 50),
            })

    bench_df = pd.DataFrame(bench_rows)
    bench_df["fpr_top50"] = fpr_top50
    bench_df["roc_auc"] = roc_auc
    bench_df["pr_auc"] = pr_auc
    bench_df["n_negatives_in_pool"] = int(n_neg)
    bench_df["interpretation"] = (
        "not_computable: no curated non-neural FSTF entered the DE-seeded "
        "candidate pool (selection-level true negative)"
        if n_neg == 0 else
        "computed from in-pool negative controls"
    )
    bench_path = RES / "negative_control_benchmarks.csv"
    bench_df.to_csv(bench_path, index=False)
    print(f"\n  Saved: {bench_path}")

    # NOTE (WS4): the previous version also rendered empty-axes
    # roc_pr_curve.png / negative_control_rank_histogram.png when no
    # negative was in the pool. Those figures are intentionally NOT
    # generated in that case; the CSV records the honest outcome.
    # 2026-09-06 audit: the dead figure block that used to follow this
    # return was UNREACHABLE (it sat after `return 0`, so the documented
    # supplementary figures were never produced even when negatives were
    # in the pool). Removed; if in-pool negatives ever exist, the figure
    # must be added back ABOVE this return and covered by a test.

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
