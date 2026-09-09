"""Evidence stream coverage across all TF candidates."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt, numpy as np

def build():
    df = load_all()
    total = len(df)
    streams = [s for s in STREAM_COLS if s in df.columns]
    counts = {s: ((df[s]>0) & df[s].notna()).sum() for s in streams}
    fracs = {s: c/total for s,c in counts.items()}
    vals = [fracs[s]*100 for s in streams]
    colors = [STREAM_C[s] for s in streams]
    labels = [STREAM_L[s] for s in streams]

    fig, ax = plt.subplots(figsize=(W_15COL, 2.8))
    y = np.arange(len(streams))[::-1]
    bars = ax.barh(y, vals, color=colors, edgecolor="none", height=0.62)
    for bar, v, c in zip(bars, vals, counts.values()):
        ax.text(bar.get_width() + 1.2, bar.get_y() + bar.get_height() / 2,
                f"{v:.1f}%", va="center", fontsize=6.2, color="#333333")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.8)
    ax.set_xlabel("Coverage across all transcription factors (%)", fontsize=7.0)
    ax.set_ylabel("Evidence stream", fontsize=7.0)
    ax.set_xlim(0, 110)
    ax.set_title(f"Evidence Stream Completeness Across Transcriptome (N = {total:,})", fontsize=8.0, pad=6)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save(fig, "01_stream_coverage_all")

if __name__=="__main__": build()

