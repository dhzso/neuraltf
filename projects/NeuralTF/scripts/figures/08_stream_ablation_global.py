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
    fig, ax = plt.subplots(figsize=(6, 4))
    y = np.arange(len(stats))
    ax.barh(y, stats["median_abs"], color=[STREAM_C[s] for s in stats["stream"]], height=0.55, edgecolor="white", lw=0.3)
    for i, (_,r) in enumerate(stats.iterrows()):
        ax.text(r["median_abs"]+0.05, i, f'ΔTop10={r["n_displaced"]}', fontsize=7, va="center",
                color=C_HL if r["n_displaced"]>0 else "#999")
    ax.set_yticks(y); ax.set_yticklabels([STREAM_L[s] for s in stats["stream"]], fontsize=8)
    ax.set_xlabel(f"Median |rank change| ({len(neural)} candidates)"); ax.set_ylabel("Evidence stream removed")
    ax.set_title("Removing expression or reproducibility causes largest rank shifts",
                 fontweight="bold", pad=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=STREAM_C[s], label=STREAM_L[s]) for s in stats["stream"]],
              loc="lower right", fontsize=6, frameon=True, title="Removed stream", title_fontsize=7)
    fig.tight_layout(); save(fig, "08_stream_ablation_global")

if __name__=="__main__": build()
