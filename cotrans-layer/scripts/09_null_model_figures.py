"""
09_null_model_figures.py
generates publication figures for the null model analysis.
figure 6: real vs null minimum gaps (scatter, 2-panel)
figure 7: p-value distributions (histogram, 2-panel)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# -- paths --
results_dir = Path("/storage/kiran-stuff/codon_project/cotrans-layer/results/summary_104")
fig_dir = Path("/storage/kiran-stuff/codon_project/cotrans-layer/results/figures_104")
fig_dir.mkdir(parents=True, exist_ok=True)

# -- style --
sns.set_style("white")
plt.rcParams.update({
    "font.size": 12,
    "figure.facecolor": "white",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# -- load data --
df = pd.read_csv(results_dir / "null_model_results.csv")
n = len(df)
print(f"loaded {n} proteins")

# ============================================================
# figure 6: real vs null minimum gaps
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

for ax, label, null_med_col, pval_col, fisher_p, frac_above in [
    (axes[0], "Null A (random boundaries)",
     "null_a_median_min_gap", "null_a_pvalue_min_gap",
     "1.41e-26", "100/104"),
    (axes[1], "Null B (shuffled sizes)",
     "null_b_median_min_gap", "null_b_pvalue_min_gap",
     "1.0", "13/104"),
]:
    x = df[null_med_col]
    y = df["real_min_gap"]
    sig = df[pval_col] < 0.05

    # non-significant points (gray)
    ax.scatter(x[~sig], y[~sig], c="gray", alpha=0.5, s=30, label="p ≥ 0.05", zorder=2)
    # significant points (firebrick)
    ax.scatter(x[sig], y[sig], c="firebrick", alpha=0.7, s=30, label="p < 0.05", zorder=3)

    # diagonal y = x
    lo = 0
    hi = max(x.max(), y.max()) * 1.1
    ax.plot([lo, hi], [lo, hi], "k--", lw=1, alpha=0.5, zorder=1)

    ax.set_xlabel("null median minimum gap (residues)")
    ax.set_ylabel("real minimum gap (residues)")
    ax.set_title(label)
    ax.legend(loc="upper left", fontsize=10, frameon=False)

    # annotation
    n_sig = sig.sum()
    txt = f"Fisher combined p = {fisher_p}\nabove diagonal: {frac_above}"
    ax.text(0.97, 0.03, txt, transform=ax.transAxes,
            ha="right", va="bottom", fontsize=9,
            bbox=dict(facecolor="white", edgecolor="gray", alpha=0.8, boxstyle="round,pad=0.3"))

fig.suptitle("Real vs null minimum L3 emergence gap", fontsize=14, y=1.02)
sns.despine(fig=fig)
plt.tight_layout()

for ext in ["png", "pdf"]:
    fig.savefig(fig_dir / f"null_model_scatter.{ext}", dpi=300, bbox_inches="tight")
    print(f"saved null_model_scatter.{ext}")
plt.close(fig)

# ============================================================
# figure 7: p-value distributions
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

bins = np.linspace(0, 1, 21)  # 20 bins from 0 to 1

for ax, pval_col, label, n_sig_expected in [
    (axes[0], "null_a_pvalue_min_gap", "Null A (random boundaries)", 28),
    (axes[1], "null_b_pvalue_min_gap", "Null B (shuffled sizes)", 0),
]:
    pvals = df[pval_col]
    n_sig = (pvals < 0.05).sum()

    ax.hist(pvals, bins=bins, color="steelblue", edgecolor="white", alpha=0.8)
    ax.axvline(0.05, color="firebrick", ls="--", lw=1.5, label="p = 0.05")
    ax.set_xlabel("p-value (min gap)")
    ax.set_ylabel("count")
    ax.set_title(label)

    # annotation with actual count from data
    ax.text(0.97, 0.95, f"significant: {n_sig}/{n}",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(facecolor="white", edgecolor="gray", alpha=0.8, boxstyle="round,pad=0.3"))

    ax.legend(loc="center right", fontsize=10, frameon=False)

fig.suptitle("P-value distribution: real minimum gap vs null", fontsize=14, y=1.02)
sns.despine(fig=fig)
plt.tight_layout()

for ext in ["png", "pdf"]:
    fig.savefig(fig_dir / f"null_model_pvalues.{ext}", dpi=300, bbox_inches="tight")
    print(f"saved null_model_pvalues.{ext}")
plt.close(fig)

print("done.")
