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

    fig, ax = plt.subplots(figsize=(6, 4))
    y = np.arange(len(streams))[::-1]
    bars = ax.barh(y, vals, color=colors, edgecolor="white", linewidth=0.5, height=0.6)
    for bar, v, c in zip(bars, vals, counts.values()):
        ax.text(bar.get_width()+0.5, bar.get_y()+bar.get_height()/2,
                f" {c}/{total} ({v:.0f}%)", va="center", fontsize=8)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(f"% of TF candidates (n={total}) with evidence"); ax.set_ylabel("Evidence stream")
    ax.set_xlim(0, 105)
    ax.set_title("Evidence stream coverage across planarian TF candidates",
                 fontweight="bold", pad=8)

    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=STREAM_C[s], label=STREAM_L[s]) for s in streams],
              loc="lower right", fontsize=6, frameon=True, title="Stream", title_fontsize=7)
    fig.tight_layout(); save(fig, "01_stream_coverage_all")

if __name__=="__main__": build()
