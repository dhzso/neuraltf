#!/usr/bin/env python
"""Bootstrap confidence intervals for integrated evidence scores.

Method (WS3): the correct unit of uncertainty is the WEIGHT VECTOR, not
the candidate rows. Resampling candidate rows mixes genes and produces a
mixture mean (~the cohort average) for every gene. Instead we consume
the per-candidate draw-score matrices written by the Dirichlet scripts
(1000 weight draws, seed 2024) and report per-gene 2.5/97.5 percentiles
of the renormalized weighted score under weight uncertainty. A parametric
stream-bootstrap fallback (resampling the observed streams with noise)
covers runs without draw matrices.

Outputs:
  results/bootstrap_scores_ci.csv  — per-gene mean + 95% CI
  figures/25_bootstrap_ci.png is produced by the numbered figure script.

Usage:
    python scripts/stats/bootstrap_confidence.py [--n-boot 1000] [--seed 42]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RUN_DIR = REPO / "projects" / "NeuralTF" / "runs" / "pipeline_run"
RESULTS_DIR = REPO / "projects" / "NeuralTF" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

STREAMS = [
    "expression", "specificity", "reproducibility", "rnai",
    "correlation", "neural_enriched", "neural_specificity",
    "perez_lineage", "perez_influence",
]
# Must match bioforge.evidence.scoring.DEFAULT_WEIGHTS exactly.
W_DEFAULT = np.array([0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])


def load_score_matrix() -> pd.DataFrame:
    candidates_path = RUN_DIR / "rank.csv"
    if not candidates_path.exists():
        raise FileNotFoundError(
            f"No candidate score file found at {candidates_path}; run the pipeline first"
        )
    return pd.read_csv(candidates_path).drop_duplicates(subset="gene_id", keep="first")


def dirichlet_draw_cis(df: pd.DataFrame) -> pd.DataFrame | None:
    """Per-gene 2.5/97.5 percentiles from the persisted Dirichlet draw
    matrices (the same weight uncertainty model the robustness analysis
    uses; no new randomness needed)."""
    out_rows = []
    for method, path in (
        ("centered", RESULTS_DIR / "dirichlet_centered_draw_scores.csv"),
        ("uniform", RESULTS_DIR / "dirichlet_uniform_draw_scores.csv"),
    ):
        if not path.exists():
            continue
        draws = pd.read_csv(path)
        gene_col = "gene_id" if "gene_id" in draws.columns else draws.columns[0]
        draw_cols = [c for c in draws.columns if c.startswith("draw_")]
        if not draw_cols:
            continue
        mat = draws[draw_cols].to_numpy(dtype=float)
        lo, hi = np.percentile(mat, [2.5, 97.5], axis=1)
        mean = mat.mean(axis=1)
        for i, gid in enumerate(draws[gene_col].astype(str)):
            out_rows.append({
                "gene_id": gid,
                f"{method}_mean": mean[i],
                f"{method}_ci_95_lo": lo[i],
                f"{method}_ci_95_hi": hi[i],
            })
    if not out_rows:
        return None
    base = pd.DataFrame(out_rows)
    # merge centered + uniform side by side on gene_id
    out = df[["gene_id", "gene_name"]].merge(
        base[base.columns[:4]], on="gene_id", how="left"
    ) if "centered_mean" in base.columns else None
    if out is None:
        return None
    uni = base[["gene_id", "uniform_mean", "uniform_ci_95_lo", "uniform_ci_95_hi"]] \
        if "uniform_mean" in base.columns else None
    if uni is not None:
        out = out.merge(uni, on="gene_id", how="left")
    return out


def weight_bootstrap_cis(df: pd.DataFrame, n_boot: int, seed: int) -> pd.DataFrame:
    """Fallback: parametric bootstrap over the weight vector.

    Draws Dirichlet(k=40 * W) weights (the centered model), recomputes the
    renormalized weighted score per gene per draw, and takes percentiles.
    This resamples the WEIGHTS (the modeling choice), never the genes.
    """
    stream_cols = [c for c in STREAMS if c in df.columns]
    S = df[stream_cols].to_numpy(dtype=float)
    W = W_DEFAULT[:len(stream_cols)]
    W = W / W.sum()

    rng = np.random.default_rng(seed)
    draws = rng.dirichlet(40.0 * W, size=n_boot)  # (n_boot, n_streams)
    valid = ~np.isnan(S)
    S_filled = np.nan_to_num(S, nan=0.0)
    num = S_filled @ draws.T                    # (n_genes, n_boot)
    den = valid.astype(float) @ draws.T
    den = np.where(den > 0, den, 1.0)
    scores = num / den

    lo, hi = np.percentile(scores, [2.5, 97.5], axis=1)
    out = df[["gene_id"] + (["gene_name"] if "gene_name" in df.columns else [])].copy()
    out["bootstrap_mean"] = scores.mean(axis=1)
    out["mean_score"] = out["bootstrap_mean"]
    out["ci_95_lo"] = lo
    out["ci_low"] = lo
    out["ci_95_hi"] = hi
    out["ci_high"] = hi
    out["ci_width"] = hi - lo
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Bootstrap 95% CI for integrated scores (weight-vector uncertainty)"
    )
    parser.add_argument("--n-boot", type=int, default=1000,
                        help="Number of weight draws (fallback mode)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (fallback mode)")
    args = parser.parse_args()

    print(f"=== Bootstrap Confidence Intervals (weight-vector model, n={args.n_boot}) ===")

    df = load_score_matrix()
    print(f"Candidates: {len(df)} (unique gene_id)")

    out = dirichlet_draw_cis(df)
    if out is not None:
        print("Using persisted Dirichlet draw matrices (centered k=40 / uniform a=1).")
    else:
        print("Draw matrices not found; running parametric weight bootstrap "
              "(Dirichlet k=40 centered at default weights).")
        out = weight_bootstrap_cis(df, args.n_boot, args.seed)

    # attach the observed integrated score for reference ordering
    if "integrated_score" in df.columns:
        out = out.merge(df[["gene_id", "integrated_score"]], on="gene_id", how="left")
        out = out.sort_values("integrated_score", ascending=False)
    out = out.reset_index(drop=True)

    out_path = RESULTS_DIR / "bootstrap_scores_ci.csv"
    out.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")

    score_col = "centered_mean" if "centered_mean" in out.columns else "bootstrap_mean"
    ci_lo_col = ("centered_ci_95_lo" if "centered_ci_95_lo" in out.columns
                 else "ci_95_lo")
    ci_hi_col = ("centered_ci_95_hi" if "centered_ci_95_hi" in out.columns
                 else "ci_95_hi")
    print(f"\nTop-10 by integrated score (mean and 95% CI under weight uncertainty):")
    for _, row in out.head(10).iterrows():
        label = row.get("gene_name", row["gene_id"])
        print(f"  {str(label):>14}  mean={row[score_col]:.4f}  "
              f"95% CI=[{row[ci_lo_col]:.4f}, {row[ci_hi_col]:.4f}]")

    return 0


if __name__ == "__main__":
    sys.exit(main())
