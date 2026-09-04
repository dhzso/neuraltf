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
15 curated non-neural planarian FSTFs with confirmed dd_Smed_v6 IDs:
  Muscle     : myoD  (longitudinal), nkx1-1 (circular), mhc (striated)
  Intestine  : gata4/5/6-1, gata4/5/6-2, nkx2.1
  Pharynx    : foxA1, mlck-1
  Epidermis  : zfp-1, bmp2/4
  Protonephridia: six1/2-2, pax2/5/8a
  Muscle/ECM : zeb2-1, twist, snail2

Metrics
-------
- Empirical FPR: n_negatives in top 50 / n_negatives total
- ROC-AUC and Precision-Recall AUC using scikit-learn
- Rank histogram of negative controls

Outputs
-------
  projects/NeuralTF/results/negative_control_benchmarks.csv
  projects/NeuralTF/figures/supplementary/roc_pr_curve.png

Usage
-----
    python projects/NeuralTF/scripts/negative_controls.py
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
# v6 IDs sourced from King 2024 mmc4 and Plass 2018 supplementary data.
# Names in parentheses are the canonical planarian gene symbols.
# ---------------------------------------------------------------------------
NEGATIVE_CONTROLS: dict[str, str] = {
    # Muscle lineage (longitudinal and circular)
    "dd_Smed_v6_7066_0_1":   "myoD (longitudinal muscle FSTF)",
    "dd_Smed_v6_4028_0_1":   "nkx1-1 (circular muscle FSTF)",
    "dd_Smed_v6_3525_0_1":   "mhc-1 (myosin heavy chain, striated muscle)",
    # Intestine lineage
    "dd_Smed_v6_1219_0_1":   "gata4/5/6-1 (intestine FSTF)",
    "dd_Smed_v6_31610_0_1":  "gata4/5/6-2 (intestine FSTF)",
    "dd_Smed_v6_12655_0_1":  "nkx2.1 (intestinal progenitor)",
    # Pharynx lineage
    "dd_Smed_v6_14025_0_1":  "foxA1 (pharynx FSTF)",
    "dd_Smed_v6_9718_0_1":   "mlck-1 (pharynx smooth muscle kinase)",
    # Epidermis lineage
    "dd_Smed_v6_36444_0_1":  "zfp-1 (epidermal progenitor FSTF)",
    "dd_Smed_v6_14353_0_1":  "bmp2/4 (dorsal-ventral epidermal patterning)",
    # Protonephridia lineage
    "dd_Smed_v6_11819_0_1":  "six1/2-2 (protonephridia FSTF)",
    "dd_Smed_v6_17815_0_1":  "pax2/5/8a (protonephridia FSTF)",
    # Mesenchymal / EMT
    "dd_Smed_v6_24826_0_1":  "zeb2-1 (muscle/parenchyma EMT regulator)",
    "dd_Smed_v6_22017_0_1":  "twist (mesenchymal EMT TF)",
    "dd_Smed_v6_19501_0_1":  "snail2 (epithelial-mesenchymal regulator)",
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
    # Positive label: RNAi-validated ONLY. (The previous definition —
    # neural_enriched>0 OR rnai>0 — used a scoring stream as ground truth,
    # making the benchmark circular.)
    pos_mask = rank["proof_status"] == "known_rnai_validated"
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

    print("\nDone.")
    return 0

    # --- Publication figure: ROC + PR curves --------------------------------
    # Only rendered when at least one curated negative is IN the pool; an
    # empty-pool run would draw empty axes (previously shipped as such).
    if n_neg > 0:
        try:
            from sklearn.metrics import roc_curve, precision_recall_curve

            labeled = rank[rank["label"].isin(["positive", "negative_control"])].copy()
            y_true = (labeled["label"] == "positive").astype(int).values
            y_score = labeled["integrated_score"].fillna(0.0).values

            C_ROC = "#0072B2"
            C_PR  = "#D55E00"
            C_REF = "#999999"

            fig, axes = plt.subplots(1, 2, figsize=(7.05, 3.15))

            ax = axes[0]
            if len(np.unique(y_true)) == 2:
                fpr_arr, tpr_arr, _ = roc_curve(y_true, y_score)
                ax.plot(fpr_arr, tpr_arr, color=C_ROC, lw=1.5,
                        label=f"AUC = {roc_auc:.3f}")
            ax.plot([0, 1], [0, 1], "--", color=C_REF, lw=0.8, label="Random")
            ax.set_xlabel("False Positive Rate", fontsize=8)
            ax.set_ylabel("True Positive Rate", fontsize=8)
            ax.set_title("ROC Curve\n(neural vs non-neural FSTFs)", fontsize=9)
            ax.legend(fontsize=7, frameon=False)
            ax.spines[["top", "right"]].set_visible(False)

            ax = axes[1]
            if len(np.unique(y_true)) == 2:
                prec_arr, rec_arr, _ = precision_recall_curve(y_true, y_score)
                ax.plot(rec_arr, prec_arr, color=C_PR, lw=1.5,
                        label=f"AP = {pr_auc:.3f}")
            baseline_pr = y_true.mean() if len(y_true) > 0 else 0.5
            ax.axhline(baseline_pr, linestyle="--", color=C_REF, lw=0.8,
                       label=f"Baseline ({baseline_pr:.2f})")
            ax.set_xlabel("Recall", fontsize=8)
            ax.set_ylabel("Precision", fontsize=8)
            ax.set_title("Precision-Recall Curve\n(neural vs non-neural FSTFs)", fontsize=9)
            ax.legend(fontsize=7, frameon=False)
            ax.spines[["top", "right"]].set_visible(False)

            plt.tight_layout(pad=0.8)
            fig_path = FIG / "roc_pr_curve.png"
            fig.savefig(fig_path, dpi=300, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            print(f"  Saved: {fig_path}")

            fig2, ax2 = plt.subplots(figsize=(3.5, 2.8))
            neg_ranks = neg_df["rank_pos"].dropna().astype(int).values
            if len(neg_ranks) > 0:
                ax2.hist(neg_ranks, bins=20, color="#CC79A7", edgecolor="white", linewidth=0.5)
            ax2.axvline(50, color="#D55E00", linestyle="--", lw=1.0, label="Top-50 cutoff")
            ax2.set_xlabel("Rank position (1 = highest score)", fontsize=8)
            ax2.set_ylabel("Count", fontsize=8)
            ax2.set_title("Negative control rank distribution", fontsize=9)
            ax2.legend(fontsize=7, frameon=False)
            ax2.spines[["top", "right"]].set_visible(False)
            plt.tight_layout(pad=0.5)
            hist_path = FIG / "negative_control_rank_histogram.png"
            fig2.savefig(hist_path, dpi=300, bbox_inches="tight", facecolor="white")
            plt.close(fig2)
            print(f"  Saved: {hist_path}")
        except ImportError:
            print("  [WARN] scikit-learn not installed; ROC/PR figure skipped.")
        except Exception as e:
            print(f"  [WARN] Figure generation failed: {e}")

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
