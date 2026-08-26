#!/usr/bin/env python
"""Regenerate every NeuralTF artifact from the raw downloads, end to end.

One command after you have placed the raw source files in their documented
locations (see README.md "Datasets" / datasets/MANIFEST.md)::

    python scripts/generate_all.py            # everything, incl. PlanMine (network)
    python scripts/generate_all.py --skip-planmine   # offline (keeps existing parquet)

Steps run in dependency order, each only if its inputs are present AND its
output does not already exist (re-run with ``--force`` to regenerate
existing outputs); a matching summary is printed at the end and missing-input
steps are skipped, never aborted:

  1. convert_fincher.py              -> datasets/processed/fincher_subsample.h5ad
  2. consolidate_plass.py            -> datasets/processed/plass_v6.h5ad
  3. build_bridge.py                 -> projects/NeuralTF/data/bridge.csv
  4. build_king_atlas.py             -> projects/NeuralTF/data/king_atlas.tsv
  5. pipeline (run.py)               -> projects/NeuralTF/runs/pipeline_run/*
6. query_planmine.py                  -> datasets/processed/planmine_annotations.parquet (network)
  7. prioritize_neural_tfs.py           -> projects/NeuralTF/results/top10_*.csv + report.md
  8. make_supp_go_figures.py            -> projects/NeuralTF/figures/supplementary/*
  9. dirichlet_prioritize.py            -> projects/NeuralTF/results/dirichlet_*.csv|md
  10. dirichlet_uniform.py              -> projects/NeuralTF/results/dirichlet_uniform_*.csv|txt
  11. dirichlet_uniform_all249.py       -> projects/NeuralTF/results/dirichlet_uniform_all249_*.csv|txt
  12. export_fstf_ranked.py             -> projects/NeuralTF/results/fstf_ranked_*.csv
  13. generate_publication_figures.py   -> projects/NeuralTF/figures/01-17_*.png (17 figs)
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def _ready(root: Path, step: str) -> str | None:
    """Return None if the step's raw inputs are present, else a reason."""
    raw = root / "datasets" / "raw"
    king = raw / "Supplementary_Data_ King_2024"
    proc = root / "datasets" / "processed"
    data = root / "projects" / "NeuralTF" / "data"
    run = root / "projects" / "NeuralTF" / "runs" / "pipeline_run"
    rules = {
        "Fincher h5ad": (
            (raw / "GSE111764_GEO_Fincher_atlas"
             / "GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz").exists(),
            "GSE111764 DGE .txt.gz missing (README 'Datasets' -> Fincher)"),
        "Plass h5ad": (
            any(p.suffix == ".tar" and ("RAW" in p.name.upper() or "GSE103633" in p.name.upper())
                for p in raw.rglob("*")) if raw.exists() else False,
            "GSE103633_RAW.tar missing (README 'Datasets' -> Plass)"),
        "Bridge CSV": (
            (raw / "smed_20140614.mapping.rosettastone.2020.txt").exists()
            and any(king.glob("*mmc4*.xlsx")) if king.exists() else False,
            "Rosetta Stone txt or King mmc4 xlsx missing"),
        "King atlas TSV": (
            any(king.glob("*mmc4*.xlsx")) and any(king.glob("*mmc7*.xlsx"))
            if king.exists() else False,
            "King mmc4/mmc7 xlsx missing"),
        "Pipeline run": (
            (proc / "fincher_subsample.h5ad").exists()
            and (proc / "plass_v6.h5ad").exists()
            and any(king.glob("*mmc4*.xlsx"))
            and any(king.glob("*mmc5*.xlsx"))
            and any(king.glob("*mmc6*.xlsx")) if king.exists() else False,
            "h5ads or King mmc4-6 xlsx missing"),
        "PlanMine parquet": (
            (run / "rank_neural.csv").exists(),
            "rank_neural.csv missing (run the pipeline first); network step"),
        "Prioritization": (
            (run / "rank_neural.csv").exists()
            and (proc / "planmine_annotations.parquet").exists()
            and any(king.glob("*mmc4*.xlsx"))
and any(king.glob("*mmc5*.xlsx")) if king.exists() else False,
            "rank_neural.csv / PlanMine parquet / bridge.csv / King mmc4-5 xlsx missing"),
        "GO supp figures": (
            (root / "projects" / "NeuralTF" / "results"
             / "top10_neural_tfs_prioritized.csv").exists()
            and (raw / "go.obo").exists(),
            "top10 shortlist (prioritization) or go.obo missing"),
        "Dirichlet-centered": (
            (run / "rank_neural.csv").exists()
            and (proc / "planmine_annotations.parquet").exists()
            and (data / "bridge.csv").exists()
            and any(king.glob("*mmc4*.xlsx"))
            and any(king.glob("*mmc5*.xlsx")) if king.exists() else False,
            "rank_neural.csv / PlanMine parquet / bridge.csv / King mmc4-5 xlsx missing"),
        "Dirichlet-uniform": (
            (run / "rank_neural.csv").exists()
            and (proc / "planmine_annotations.parquet").exists()
            and (data / "bridge.csv").exists()
            and any(king.glob("*mmc4*.xlsx"))
            and any(king.glob("*mmc5*.xlsx")) if king.exists() else False,
            "rank_neural.csv / PlanMine parquet / bridge.csv / King mmc4-5 xlsx missing"),
        "Dirichlet-all249": (
            (run / "rank.csv").exists()
            and (run / "rank_neural.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "dirichlet_uniform_top10.csv").exists(),
            "rank.csv / rank_neural.csv / dirichlet_uniform_top10.csv missing"),
        "Export ranked FSTF": (
            (run / "rank.csv").exists(),
            "rank.csv missing (run the pipeline first)"),
        "Visualize fixed-weight": (
            (run / "rank.csv").exists()
            and (run / "rank_neural.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "top10_neural_tfs_prioritized.csv").exists(),
            "rank.csv / rank_neural.csv / top10_neural_tfs_prioritized.csv missing"),
        "Visualize Dirichlet-centered": (
            (root / "projects" / "NeuralTF" / "results"
             / "dirichlet_top10_prioritized.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "top10_neural_tfs_prioritized.csv").exists(),
            "dirichlet_top10_prioritized.csv / top10_neural_tfs_prioritized.csv missing"),
        "Visualize Dirichlet-uniform": (
            (root / "projects" / "NeuralTF" / "results"
             / "dirichlet_uniform_top10.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "dirichlet_top10_prioritized.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "top10_neural_tfs_prioritized.csv").exists(),
            "dirichlet_uniform_top10.csv / dirichlet_top10_prioritized.csv / top10_neural_tfs_prioritized.csv missing"),
        "Visualize method comparison": (
            (root / "projects" / "NeuralTF" / "results"
             / "dirichlet_uniform_full_rank.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "dirichlet_top10_prioritized.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "dirichlet_uniform_all249_full_rank.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "top10_neural_tfs_prioritized.csv").exists(),
            "dirichlet_uniform_full_rank.csv / dirichlet_top10_prioritized.csv / "
            "dirichlet_uniform_all249_full_rank.csv / top10_neural_tfs_prioritized.csv missing"),
        "Generate supplementary tables": (
            (root / "projects" / "NeuralTF" / "results"
             / "top10_neural_tfs_prioritized.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "dirichlet_top10_prioritized.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "dirichlet_uniform_top10.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "dirichlet_uniform_all249_full_rank.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "fstf_ranked_19_neural.csv").exists(),
            "prioritization CSVs / FSTF CSVs missing"),
        "Generate publication figures": (
            (run / "rank.csv").exists()
            and (run / "rank_neural.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "top10_neural_tfs_prioritized.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "dirichlet_top10_prioritized.csv").exists()
            and (root / "projects" / "NeuralTF" / "results"
                 / "dirichlet_uniform_top10.csv").exists(),
            "rank.csv / rank_neural.csv / prioritization CSVs missing"),
        "Validate with Perez ANANSE": (
            (root / "projects" / "NeuralTF" / "results"
             / "top10_neural_tfs_prioritized.csv").exists()
            and (root / "datasets" / "raw" / "Supplementary_Data_ Perez_2025"
                 / "41467_2025_65712_MOESM22_ESM.xlsx").exists(),
            "prioritized top10 or Perez MOESM22 missing"),
    }
    ok, why = rules[step]
    return None if ok else why


STEPS = [
    ("Fincher h5ad", ["scripts", "convert_fincher.py"], [],
     [["datasets", "processed", "fincher_subsample.h5ad"]]),
    ("Plass h5ad", ["scripts", "consolidate_plass.py"], [],
     [["datasets", "processed", "plass_v6.h5ad"]]),
    ("Bridge CSV", ["scripts", "build_bridge.py"], [],
     [["projects", "NeuralTF", "data", "bridge.csv"]]),
    ("King atlas TSV", ["scripts", "build_king_atlas.py"], [],
     [["projects", "NeuralTF", "data", "king_atlas.tsv"]]),
    ("Pipeline run", ["scripts", "run.py"], [],
     [["projects", "NeuralTF", "runs", "pipeline_run", "rank.csv"],
      ["projects", "NeuralTF", "runs", "pipeline_run", "rank_neural.csv"],
      ["projects", "NeuralTF", "runs", "pipeline_run", "pipeline_results.json"],
      ["projects", "NeuralTF", "runs", "pipeline_run", "evidence_cards.md"]]),
    ("PlanMine parquet", ["scripts", "query_planmine.py"], ["--repo", str(REPO)],
     [["datasets", "processed", "planmine_annotations.parquet"]]),
    ("Prioritization", ["scripts", "prioritize_neural_tfs.py"], ["--repo", str(REPO)],
     [["projects", "NeuralTF", "results", "top10_neural_tfs_prioritized.csv"],
      ["projects", "NeuralTF", "results", "candidate_summary_report.md"]]),
    ("GO supp figures", ["projects/NeuralTF/scripts/make_supp_go_figures.py"], [],
     [["projects", "NeuralTF", "figures", "supplementary", "fig_s1_go_gene_term_map.png"],
      ["projects", "NeuralTF", "figures", "supplementary", "fig_s4_go_neural_focus.png"],
      ["projects", "NeuralTF", "figures", "supplementary",
       "go_gene_term_matrix_reduced.csv"],
      ["projects", "NeuralTF", "figures", "go_term_reference.csv"]]),
    ("Dirichlet-centered", ["projects/NeuralTF/scripts/dirichlet_prioritize.py"], [],
     [["projects", "NeuralTF", "results", "dirichlet_top10_prioritized.csv"],
      ["projects", "NeuralTF", "results", "dirichlet_overall_top10.csv"],
      ["projects", "NeuralTF", "results", "dirichlet_overall_top10_byscore.csv"],
      ["projects", "NeuralTF", "results", "dirichlet_candidate_summary_report.md"]]),
    ("Dirichlet-uniform", ["projects/NeuralTF/scripts/dirichlet_uniform.py"], [],
     [["projects", "NeuralTF", "results", "dirichlet_uniform_top10.csv"],
      ["projects", "NeuralTF", "results", "dirichlet_uniform_overall_top10.csv"],
      ["projects", "NeuralTF", "results", "dirichlet_uniform_full_rank.csv"],
      ["projects", "NeuralTF", "results", "dirichlet_uniform_summary.txt"]]),
    ("Dirichlet-all249",
     ["projects/NeuralTF/scripts/dirichlet_uniform_all249.py"], [],
     [["projects", "NeuralTF", "results", "dirichlet_uniform_all249_top10.csv"],
      ["projects", "NeuralTF", "results", "dirichlet_uniform_all249_overall_top10.csv"],
      ["projects", "NeuralTF", "results", "dirichlet_uniform_all249_full_rank.csv"],
      ["projects", "NeuralTF", "results", "dirichlet_uniform_all249_summary.txt"]]),
    ("Export ranked FSTF",
     ["projects/NeuralTF/scripts/export_fstf_ranked.py"], [],
     [["projects", "NeuralTF", "results", "fstf_ranked_19_neural.csv"],
      ["projects", "NeuralTF", "results", "fstf_ranked_43_all.csv"],
      ["projects", "NeuralTF", "results", "fstf_ranked_74_catalog.csv"]]),
    ("Generate publication figures",
     ["projects/NeuralTF/scripts/generate_publication_figures.py"], [],
     [["projects", "NeuralTF", "figures", "01_stream_coverage_249.png"],
      ["projects", "NeuralTF", "figures", "02_integrated_vs_composite.png"],
      ["projects", "NeuralTF", "figures", "03_score_distribution_249_vs_99.png"],
      ["projects", "NeuralTF", "figures", "04_evidence_heatmap_99.png"],
      ["projects", "NeuralTF", "figures", "05_top10_candidate_atlas.png"],
      ["projects", "NeuralTF", "figures", "06_weight_sensitivity_ranks.png"],
      ["projects", "NeuralTF", "figures", "07_weight_sensitivity_ptop10.png"],
      ["projects", "NeuralTF", "figures", "08_stream_ablation_global.png"],
      ["projects", "NeuralTF", "figures", "09_stream_ablation_candidate.png"],
      ["projects", "NeuralTF", "figures", "10_centered_top10_scores.png"],
      ["projects", "NeuralTF", "figures", "11_centered_scatter_99.png"],
      ["projects", "NeuralTF", "figures", "12_uniform_top10_scores.png"],
      ["projects", "NeuralTF", "figures", "13_uniform_scatter_99.png"],
      ["projects", "NeuralTF", "figures", "14_uniform_99vs249_rankrank.png"],
      ["projects", "NeuralTF", "figures", "15_method_bumpchart.png"],
      ["projects", "NeuralTF", "figures", "16_method_score_density.png"],
      ["projects", "NeuralTF", "figures", "17_method_rank_correlation.png"],
      ["projects", "NeuralTF", "figures", "18_composite_bonus_waterfall.png"],
      ["projects", "NeuralTF", "figures", "19_method_consensus.png"],
      ["projects", "NeuralTF", "figures", "20_stream_correlation.png"],
      ["projects", "NeuralTF", "figures", "21_centered_vs_uniform_scatter.png"]]),
    ("Validate with Perez ANANSE",
     ["projects/NeuralTF/scripts/validate_with_perez.py"], [],
     []),  # validation script, no persistent outputs
    ("Generate supplementary tables",
     ["projects/NeuralTF/scripts/create_supplementary_tables.py"], [],
     [["projects", "NeuralTF", "results", "supplementary_table_S1_method_comparison.csv"],
      ["projects", "NeuralTF", "results", "supplementary_table_S2_fixed_all249.csv"],
      ["projects", "NeuralTF", "results", "supplementary_table_S3_centered_all99.csv"],
      ["projects", "NeuralTF", "results", "supplementary_table_S4_uniform_all99.csv"],
      ["projects", "NeuralTF", "results", "supplementary_table_S5_uniform_all249.csv"],
      ["projects", "NeuralTF", "results", "supplementary_table_S6_fstf_19_neural.csv"],
      ["projects", "NeuralTF", "results", "supplementary_table_S7_fstf_43_all.csv"],
      ["projects", "NeuralTF", "results", "supplementary_table_S8_fstf_74_catalog.csv"],
      ["projects", "NeuralTF", "results", "supplementary_table_S9_99vs249.csv"]]),
]


def _has_output(root: Path, out_parts: list[list[str]]) -> bool:
    """True only when *every* expected output of a step exists and is non-empty
    (a crashed run that wrote one of several outputs must not be skipped).
    Empty output list means no outputs to check — must run the step."""
    if not out_parts:
        return False  # No outputs defined → must run
    for parts in out_parts:
        out = root.joinpath(*parts)
        if not (out.exists() and out.stat().st_size > 0):
            return False
    return True


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--force", action="store_true",
                   help="regenerate outputs that already exist (default: skip them)")
    p.add_argument("--skip-planmine", action="store_true",
                   help="skip the PlanMine network query (keeps existing parquet)")
    p.add_argument("--repo", type=Path, default=REPO, help="repo root (default: this repo)")
    args = p.parse_args()

    root = args.repo.resolve()
    print(f"== regenerate-all: {root} ==", flush=True)
    results: list[tuple[str, str]] = []
    for name, rel_script, extra, out_parts in STEPS:
        if name == "PlanMine parquet" and args.skip_planmine:
            results.append((name, "skipped (--skip-planmine)"))
            continue
        if not args.force and _has_output(root, out_parts):
            results.append((name, "skipped (output exists; --force to regenerate)"))
            continue
        why = _ready(root, name)
        if why:
            results.append((name, f"skipped - {why}"))
            continue
        script = root.joinpath(*rel_script)
        print(f"\n>>> {name}  [python {script.relative_to(root)}]", flush=True)
        try:
            res = subprocess.run(
                [sys.executable, str(script), *extra], cwd=root,
                timeout=3600 * 4,
            )
        except subprocess.TimeoutExpired:
            results.append((name, "FAILED (timeout)"))
            continue
        results.append((name, "OK" if res.returncode == 0 else f"FAILED (exit {res.returncode})"))

    print("\n=== summary ===", flush=True)
    for name, status in results:
        print(f"  {name:<20} {status}", flush=True)
    failed = [r for r in results if r[1].startswith("FAILED")]
    if failed:
        print(f"\n{len(failed)} step(s) failed - see output above.", flush=True)
    else:
        print("\nAll steps completed successfully.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())