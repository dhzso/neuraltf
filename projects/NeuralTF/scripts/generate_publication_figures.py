"""Master script — regenerate all publication figures (PNG only, 300dpi).

Usage:
    python projects/NeuralTF/scripts/generate_publication_figures.py
    python projects/NeuralTF/scripts/generate_publication_figures.py --figure 1 5 10
"""
from __future__ import annotations
import argparse, importlib.util, sys, time, traceback
from pathlib import Path

FIGURES_DIR = Path(__file__).resolve().parent / "figures"

FIGURES = {
    1:  ("01_stream_coverage_all.py",            "Evidence stream coverage (all TFs)"),
    3:  ("03_score_distribution_all_vs_neural.py", "Score distribution all vs neural"),
    4:  ("04_evidence_heatmap_neural.py",         "Evidence heatmap (neural)"),
    5:  ("05_top10_candidate_atlas.py",           "Top 10 candidate atlas"),
    6:  ("06_weight_sensitivity_ranks.py",        "Weight sensitivity rank distributions"),
    7:  ("07_weight_sensitivity_ptop10.py",       "Weight sensitivity P(Top10)"),
    8:  ("08_stream_ablation_global.py",          "Stream ablation global impact"),
    9:  ("09_stream_ablation_candidate.py",       "Stream ablation candidate sensitivity"),
    13: ("13_uniform_scatter_all.py",             "Fixed vs uniform Dirichlet (all)"),
    14: ("14_uniform_neural_vs_all_rankrank.py",  "Neural vs all rank-rank comparison"),
    15: ("15_method_bumpchart.py",                "3-method rank comparison"),
    17: ("17_method_rank_correlation.py",         "3-method rank correlation"),
    18: ("18_composite_bonus_waterfall.py",       "Composite bonus waterfall"),
    20: ("20_stream_correlation.py",              "Stream correlation matrix"),
    22: ("22_pipeline_schematic.py",              "Pipeline schematic"),
    23: ("23_roc_pr_curve.py",                    "ROC and PR curves"),
    24: ("24_negative_controls.py",               "Negative controls"),
    25: ("25_bootstrap_ci.py",                    "Bootstrap confidence intervals"),
    26: ("26_permutation_null.py",                "Permutation null distribution"),
    28: ("28_effect_sizes.py",                    "Effect sizes"),
    29: ("29_convergence_analysis.py",            "Convergence analysis"),
    30: ("30_calibration.py",                     "Calibration plot"),
    31: ("31_score_distribution_all9.py",         "Score distribution (9 streams)"),
    32: ("32_perez_influence_comparison.py",      "Perez influence comparison"),
    33: ("33_method_agreement_summary.py",        "Method agreement summary"),
    34: ("34_ananse_regulatory_network.py",       "ANANSE regulatory network & targets"),
    35: ("35_meta_analysis_concordance.py",       "Cross-atlas meta-analysis concordance"),
}

COMPOSITE_SCRIPTS = [
    ("composite_nature_figures.py", "Composite Nature manuscript figures (Figs 1-5)"),
    ("supp_go_figures.py",          "Supplementary GO figures (Figs S1-S2)"),
]

def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--figure", nargs="*", type=int)
    p.add_argument("--composite-only", action="store_true", help="Only run composite Nature figures")
    args = p.parse_args()

    t0 = time.time()
    ok, fail = [], []

    if not args.composite_only:
        nums = args.figure if args.figure else sorted(FIGURES.keys())
        for num in nums:
            if num not in FIGURES:
                print(f"  [SKIP] Figure {num}"); continue
            fname, desc = FIGURES[num]
            print(f"\n  [{num}] {desc}")
            try:
                mod = _load(FIGURES_DIR / fname, f"fig{num}")
                mod.build()
                ok.append(num)
                print(f"    OK")
            except Exception:
                traceback.print_exc()
                fail.append(num)

    if not args.figure or args.composite_only:
        for fname, desc in COMPOSITE_SCRIPTS:
            print(f"\n  [COMPOSITE/SUPP] {desc}")
            try:
                mod = _load(FIGURES_DIR / fname, fname.replace(".py", ""))
                if hasattr(mod, "build_all"):
                    mod.build_all()
                elif hasattr(mod, "build"):
                    mod.build()
                elif hasattr(mod, "fig_s1_go_landscape"):
                    mod.fig_s1_go_landscape()
                    mod.fig_s2_go_namespace_and_track()
                ok.append(fname)
                print("    OK")
            except Exception:
                traceback.print_exc()
                fail.append(fname)

    print(f"\n{'='*50}")
    print(f"  Done: {len(ok)} succeeded, {len(fail)} failed, {time.time()-t0:.1f}s")
    if fail: print(f"  Failed (non-fatal): {fail}")

if __name__ == "__main__":
    main()
