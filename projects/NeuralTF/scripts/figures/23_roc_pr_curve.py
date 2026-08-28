"""ROC and PR curves using RNAi TFs as ground truth."""
from __future__ import annotations
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from style import *
import matplotlib.pyplot as plt
import numpy as np
import json
from sklearn.metrics import roc_curve, auc, precision_recall_curve

def build():
    try:
        data_path = RES / "precision_recall.json"
        with open(data_path) as f:
            data = json.load(f)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

        # ROC curve
        fpr = np.array(data["roc"]["fpr"])
        tpr = np.array(data["roc"]["tpr"])
        roc_auc = auc(fpr, tpr)
        ax1.plot(fpr, tpr, color=C_A, lw=2, label=f"Integrated score (AUC = {roc_auc:.3f})")
        ax1.plot([0, 1], [0, 1], color="#999999", lw=1, linestyle="--", label="Random (AUC = 0.5)")
        ax1.set_xlabel("False Positive Rate")
        ax1.set_ylabel("True Positive Rate")
        ax1.set_title("ROC curve: integrated score\nvs RNAi-validated TFs", fontweight="bold", fontsize=10)
        ax1.legend(loc="lower right", fontsize=8)
        ax1.set_xlim([0, 1])
        ax1.set_ylim([0, 1.02])
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)

        # PR curve
        precision = np.array(data["pr"]["precision"])
        recall = np.array(data["pr"]["recall"])
        pr_auc = auc(recall, precision)
        ax2.plot(recall, precision, color=C_B, lw=2, label=f"Integrated score (AP = {pr_auc:.3f})")
        baseline = data["pr"].get("baseline", 0.15)
        ax2.axhline(y=baseline, color="#999999", lw=1, linestyle="--", label=f"Random baseline ({baseline:.2f})")
        ax2.set_xlabel("Recall")
        ax2.set_ylabel("Precision")
        ax2.set_title("PR curve: integrated score\nvs RNAi-validated TFs", fontweight="bold", fontsize=10)
        ax2.legend(loc="upper right", fontsize=8)
        ax2.set_xlim([0, 1])
        ax2.set_ylim([0, 1.02])
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)

        fig.tight_layout()
        save(fig, "23_roc_pr_curve")
    except FileNotFoundError as e:
        print(f"  [SKIP] {__file__}: {e}")
        return

if __name__ == "__main__":
    build()
