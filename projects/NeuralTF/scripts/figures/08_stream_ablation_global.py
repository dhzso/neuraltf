"""Stream ablation — global impact of removing each evidence stream."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np, pandas as pd

def _ablate(df, exclude_idx):
    scores = np.zeros(len(df))
    for i, (_, row) in enumerate(df.iterrows()):
        vals = np.array([row.get(s, np.nan) for s in STREAM_COLS], dtype=float)
        w = W.copy(); w[exclude_idx] = 0
        present = ~np.isnan(vals) & (vals != 0)
        if present.any():
            scores[i] = np.sum(w[present]*vals[present]) / np.sum(w[present])
    return scores

def build():
    neural = load_neural()
    base_scores = _ablate(neural, -1)  # all streams (none excluded)
    base_ranks = pd.Series(base_scores).rank(ascending=False).values

    stats = []
    for j, s in enumerate(STREAM_COLS):
        ab_scores = _ablate(neural, j)
        ab_ranks = pd.Series(ab_scores).rank(ascending=False).values
        delta = ab_ranks - base_ranks
        n_displaced = ((base_ranks<=10) & (ab_ranks>10)).sum()
        stats.append({"stream":s, "median_abs":np.median(np.abs(delta)),
                      "mean_abs":np.mean(np.abs(delta)), "n_displaced":n_displaced})

    stats = pd.DataFrame(stats).sort_values("median_abs", ascending=True)
    fig, ax = plt.subplots(figsize=(W_15COL, 2.8))
    y = np.arange(len(stats))
    ax.barh(y, stats["median_abs"], color=[STREAM_C[s] for s in stats["stream"]],
            height=0.62, edgecolor="none")
    for i, (_, r) in enumerate(stats.iterrows()):
        nd = int(r["n_displaced"])
        txt = f"{r['median_abs']:.1f} ({nd} displaced)" if nd > 0 else f"{r['median_abs']:.1f}"
        ax.text(r["median_abs"] + 0.08, i, txt, fontsize=6.0, va="center",
                color="#333333")
        
    ax.set_yticks(y)
    ax.set_yticklabels([STREAM_L[s] for s in stats["stream"]], fontsize=6.5)
    ax.set_xlabel("Median |Δrank|", fontsize=7.0)
    ax.set_ylabel("Omitted stream", fontsize=7.0)
    ax.set_xlim(0, max(stats["median_abs"]) * 1.35)
    ax.set_title("Global stream ablation impact", fontsize=8.0, pad=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save(fig, "08_stream_ablation_global")

if __name__=="__main__": build()

