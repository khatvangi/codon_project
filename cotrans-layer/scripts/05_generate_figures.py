#!/usr/bin/env python3
"""
05_generate_figures.py — publication-quality figures for cotranslational folding analysis.

generates 4 figures:
  1. emergence gap distribution (histogram)
  2. minimum kf required vs empirical range (scatter + band)
  3. emergence gap vs segment separation (scatter + trend)
  4. per-protein folding timelines (4 example proteins)

outputs:
  results/figures/emergence_gap_distribution.{png,pdf}
  results/figures/kf_required_vs_empirical.{png,pdf}
  results/figures/gap_vs_separation.{png,pdf}
  results/figures/folding_timelines.{png,pdf}

usage:
  python scripts/05_generate_figures.py
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from pathlib import Path
from scipy import stats

BASE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE))

# -- paths --
GAPS_CSV = BASE / "results" / "summary" / "emergence_gaps.csv"
KF_CSV = BASE / "results" / "summary" / "min_kf_required.csv"
SEGTYPES_CSV = BASE / "data" / "upstream" / "segment_types.csv"
FIG_DIR = BASE / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# -- global style --
plt.rcParams.update({
    "font.size": 12,
    "axes.grid": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 300,
    "font.family": "sans-serif",
})
sns.set_style("white")

# -- translation rate constant: 0.06 seconds per residue --
SECONDS_PER_RESIDUE = 0.06


def save_fig(fig, name):
    """save figure as both png and pdf."""
    fig.savefig(FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    print(f"  saved {name}.png and {name}.pdf")


def figure1_emergence_gap_distribution(gaps_df):
    """histogram of emergence gaps for all L3 interfaces."""
    print("figure 1: emergence gap distribution")

    gap_res = gaps_df["emergence_gap_residues"].values
    n = len(gap_res)
    min_gap = gap_res.min()
    median_gap = np.median(gap_res)

    fig, ax1 = plt.subplots(figsize=(8, 5))

    # histogram
    bins = np.arange(0, gap_res.max() + 10, 5)
    ax1.hist(gap_res, bins=bins, color="steelblue", edgecolor="white", linewidth=0.5)

    # vertical dashed line at minimum gap
    ax1.axvline(min_gap, color="firebrick", linestyle="--", linewidth=1.5, label=f"min = {min_gap:.0f} res")

    ax1.set_xlabel("Emergence gap (residues)")
    ax1.set_ylabel("Number of L3 interfaces")

    # secondary x-axis on top: gap in seconds
    ax2 = ax1.secondary_xaxis("top", functions=(
        lambda x: x * SECONDS_PER_RESIDUE,
        lambda x: x / SECONDS_PER_RESIDUE
    ))
    ax2.set_xlabel("Emergence gap (seconds)")

    # annotation
    min_sec = min_gap * SECONDS_PER_RESIDUE
    med_sec = median_gap * SECONDS_PER_RESIDUE
    annotation = (
        f"N = {n:,} interfaces\n"
        f"Min gap = {min_gap:.0f} residues ({min_sec:.1f} s)\n"
        f"Median = {median_gap:.0f} residues ({med_sec:.1f} s)"
    )
    ax1.text(0.97, 0.95, annotation, transform=ax1.transAxes,
             ha="right", va="top", fontsize=11,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="gray", alpha=0.9))

    sns.despine(ax=ax1)
    fig.tight_layout()
    save_fig(fig, "emergence_gap_distribution")
    plt.close(fig)


def figure2_kf_required_vs_empirical(kf_df):
    """scatter of kf_min_required vs segment size, with empirical range band."""
    print("figure 2: kf required vs empirical range")

    fig, ax = plt.subplots(figsize=(8, 5.5))

    # sort by segment size for smooth fill_between
    sizes = kf_df["seg_size"].values
    size_range = np.arange(sizes.min(), sizes.max() + 1)

    # build empirical band: for each segment size, get the empirical range
    # the empirical range depends on size via the Plaxco correlation
    # use the actual data points to define the band
    kf_sorted = kf_df.sort_values("seg_size")

    # group by seg_size and get empirical bounds (take the range that spans all entries)
    grouped = kf_sorted.groupby("seg_size").agg(
        emp_low=("empirical_kf_range_log10_low", "min"),
        emp_high=("empirical_kf_range_log10_high", "max"),
    ).reset_index()

    # fill between empirical range — use step-like appearance for clarity
    ax.fill_between(grouped["seg_size"], grouped["emp_low"], grouped["emp_high"],
                    color="lightgray", alpha=0.5, label="Empirical $k_f$ range", zorder=1)

    # plot boundaries of empirical range as lines
    ax.plot(grouped["seg_size"], grouped["emp_low"], color="gray", linewidth=0.8, linestyle="-", zorder=2)
    ax.plot(grouped["seg_size"], grouped["emp_high"], color="gray", linewidth=0.8, linestyle="-", zorder=2)

    # scatter: kf_min_required
    ax.scatter(kf_df["seg_size"], kf_df["kf_min_log10"],
               color="firebrick", s=30, alpha=0.7, edgecolors="darkred", linewidth=0.3,
               label="Min $k_f$ required", zorder=3)

    ax.set_xlabel("Foldon size (residues)")
    ax.set_ylabel(r"$\log_{10}(k_f$ / s$^{-1})$")

    # count within range
    n_within = kf_df["within_empirical_range"].sum()
    n_total = len(kf_df)
    ax.text(0.97, 0.05, f"{n_within}/{n_total} ({100*n_within/n_total:.0f}%) within empirical range",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=11,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="gray", alpha=0.9))

    ax.legend(loc="upper right", framealpha=0.9)
    sns.despine(ax=ax)
    fig.tight_layout()
    save_fig(fig, "kf_required_vs_empirical")
    plt.close(fig)


def figure3_gap_vs_separation(gaps_df):
    """scatter of emergence gap vs segment index separation."""
    print("figure 3: gap vs segment separation")

    # compute segment separation
    df = gaps_df.copy()
    df["seg_separation"] = (df["seg_j"] - df["seg_i"]).abs()

    fig, ax = plt.subplots(figsize=(8, 5.5))

    # assign color per protein
    proteins = df["protein_id"].unique()
    cmap = plt.cm.tab20
    colors = {p: cmap(i % 20) for i, p in enumerate(sorted(proteins))}

    for pid in sorted(proteins):
        sub = df[df["protein_id"] == pid]
        ax.scatter(sub["seg_separation"], sub["emergence_gap_residues"],
                   color=colors[pid], s=18, alpha=0.5, edgecolors="none")

    # trend line (linear regression)
    x = df["seg_separation"].values
    y = df["emergence_gap_residues"].values
    slope, intercept, r_val, p_val, std_err = stats.linregress(x, y)
    x_fit = np.linspace(x.min(), x.max(), 100)
    y_fit = slope * x_fit + intercept
    ax.plot(x_fit, y_fit, color="black", linewidth=1.5, linestyle="--", zorder=5)

    ax.set_xlabel("Segment index separation")
    ax.set_ylabel("Emergence gap (residues)")

    # annotate with pearson r
    p_str = f"{p_val:.2e}" if p_val > 0 else "< 1e-300"
    annotation = f"Pearson r = {r_val:.3f}\np {p_str}\nn = {len(df):,}"
    ax.text(0.03, 0.95, annotation, transform=ax.transAxes,
            ha="left", va="top", fontsize=11,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="gray", alpha=0.9))

    sns.despine(ax=ax)
    fig.tight_layout()
    save_fig(fig, "gap_vs_separation")
    plt.close(fig)


def figure4_folding_timelines(seg_df):
    """horizontal timeline diagrams for 4 example proteins."""
    print("figure 4: folding timelines")

    # tunnel length in residues (~30 AA ≈ ribosome exit tunnel)
    TUNNEL_LEN = 30

    # pick 4 example proteins: F01 (small), F38 (ubiquitin, mostly B),
    # barnase (medium, mostly A), F02 (large, 14 foldons)
    example_proteins = ["F01", "F38", "barnase", "F02"]

    fig, axes = plt.subplots(len(example_proteins), 1, figsize=(10, 8))
    fig.subplots_adjust(hspace=0.55)

    # colors for segment types
    type_colors = {"A": "steelblue", "B": "salmon"}

    for idx, pid in enumerate(example_proteins):
        ax = axes[idx]
        sub = seg_df[seg_df["protein_id"] == pid].sort_values("seg_idx")
        pdb_id = sub["pdb_id"].iloc[0]
        n_res = sub["n_res"].iloc[0]

        # y position for the foldon bar
        bar_y = 0.5
        bar_height = 0.35

        for _, row in sub.iterrows():
            seg_start = row["seg_start"]
            seg_end = row["seg_end"]
            seg_type = row["seg_type"]
            seg_idx = row["seg_idx"]

            # draw foldon bar
            color = type_colors.get(seg_type, "gray")
            ax.barh(bar_y, seg_end - seg_start, left=seg_start, height=bar_height,
                    color=color, edgecolor="white", linewidth=0.8, zorder=2)

            # label foldon with seg_idx
            mid_x = (seg_start + seg_end) / 2
            ax.text(mid_x, bar_y, str(seg_idx), ha="center", va="center",
                    fontsize=8, fontweight="bold", color="white", zorder=3)

            # emergence point: seg_end + tunnel
            emerge_x = seg_end + TUNNEL_LEN
            if emerge_x <= n_res + TUNNEL_LEN + 10:
                ax.axvline(emerge_x, color=color, linestyle=":", linewidth=0.8,
                           alpha=0.6, ymin=0.0, ymax=1.0, zorder=1)
                # small triangle marker at top
                ax.plot(emerge_x, bar_y + bar_height / 2 + 0.08, marker="v",
                        color=color, markersize=5, zorder=4)

        # axis formatting
        ax.set_xlim(-5, max(n_res + TUNNEL_LEN + 15, n_res * 1.15))
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xlabel("Residue position" if idx == len(example_proteins) - 1 else "")
        ax.set_title(f"{pid}  ({pdb_id.upper()}, {n_res} residues, {len(sub)} foldons)",
                     fontsize=11, loc="left", fontweight="bold")
        sns.despine(ax=ax, left=True)

    # add a shared legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="steelblue", edgecolor="white", label="Type A (autonomous)"),
        Patch(facecolor="salmon", edgecolor="white", label="Type B (cooperative)"),
    ]
    axes[0].legend(handles=legend_elements, loc="upper right", fontsize=9, framealpha=0.9)

    save_fig(fig, "folding_timelines")
    plt.close(fig)


def main():
    print("loading data...")
    gaps_df = pd.read_csv(GAPS_CSV)
    kf_df = pd.read_csv(KF_CSV)
    seg_df = pd.read_csv(SEGTYPES_CSV)

    print(f"  emergence gaps: {len(gaps_df)} L3 interfaces")
    print(f"  kf required: {len(kf_df)} type A foldons")
    print(f"  segment types: {len(seg_df)} segments across {seg_df['protein_id'].nunique()} proteins")
    print()

    figure1_emergence_gap_distribution(gaps_df)
    figure2_kf_required_vs_empirical(kf_df)
    figure3_gap_vs_separation(gaps_df)
    figure4_folding_timelines(seg_df)

    print()
    print(f"all figures saved to {FIG_DIR}/")


if __name__ == "__main__":
    main()
