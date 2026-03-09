#!/usr/bin/env python3
"""
07_generate_figures_104.py — publication figures for the 104-protein analysis.

generates 5 figures:
  1. emergence gap distribution — histogram with 37-protein overlay
  2. kf required vs empirical — scatter with Galpern/validation coloring
  3. gap vs segment separation — colored by protein source
  4. folding timelines — 6 example proteins (3 Galpern + 3 validation)
  5. galpern vs validation comparison — box/violin plot

outputs → results/figures_104/

usage:
  python scripts/07_generate_figures_104.py
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
GAPS_CSV = BASE / "results" / "summary_104" / "emergence_gaps.csv"
KF_CSV = BASE / "results" / "summary_104" / "min_kf_required.csv"
SEGTYPES_CSV = BASE / "data" / "extended_segment_types.csv"

# original 37-protein results for comparison
GAPS_37_CSV = BASE / "results" / "summary" / "emergence_gaps.csv"

FIG_DIR = BASE / "results" / "figures_104"
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

SECONDS_PER_RESIDUE = 0.06

# colors
GALPERN_COLOR = "steelblue"
VALIDATION_COLOR = "darkorange"
TYPE_A_COLOR = "steelblue"
TYPE_B_COLOR = "salmon"


def is_galpern(pid):
    """check if protein_id is from the Galpern set (starts with F)."""
    return pid.startswith("F")


def save_fig(fig, name):
    fig.savefig(FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    print(f"  saved {name}.png and {name}.pdf")


def figure1_emergence_gap_distribution(gaps_df):
    """histogram of emergence gaps with 37-protein reference."""
    print("figure 1: emergence gap distribution (104 proteins)")

    gap_res = gaps_df["emergence_gap_residues"].values

    # load 37-protein gaps for overlay
    try:
        gaps_37 = pd.read_csv(GAPS_37_CSV)
        gap_37_res = gaps_37["emergence_gap_residues"].values
        has_37 = True
    except FileNotFoundError:
        has_37 = False

    fig, ax1 = plt.subplots(figsize=(8, 5))

    bins = np.arange(0, gap_res.max() + 10, 5)

    # 104-protein histogram
    ax1.hist(gap_res, bins=bins, color="steelblue", edgecolor="white",
             linewidth=0.5, alpha=0.85, label=f"104 proteins (N={len(gap_res):,})")

    # 37-protein overlay (outline only)
    if has_37:
        ax1.hist(gap_37_res, bins=bins, edgecolor="firebrick",
                 linewidth=1.5, histtype="step", label=f"37 Galpern (N={len(gap_37_res):,})")

    # vertical lines
    min_gap = gap_res.min()
    ax1.axvline(min_gap, color="firebrick", linestyle="--", linewidth=1.5,
                label=f"min = {min_gap:.0f} res ({min_gap * SECONDS_PER_RESIDUE:.1f} s)")

    ax1.set_xlabel("Emergence gap (residues)")
    ax1.set_ylabel("Number of L3 interfaces")

    # secondary x-axis: seconds
    ax2 = ax1.secondary_xaxis("top", functions=(
        lambda x: x * SECONDS_PER_RESIDUE,
        lambda x: x / SECONDS_PER_RESIDUE
    ))
    ax2.set_xlabel("Emergence gap (seconds)")

    # annotation
    median_gap = np.median(gap_res)
    annotation = (
        f"N = {len(gap_res):,} interfaces\n"
        f"Min = {min_gap:.0f} res ({min_gap * SECONDS_PER_RESIDUE:.1f} s)\n"
        f"Median = {median_gap:.0f} res ({median_gap * SECONDS_PER_RESIDUE:.1f} s)\n"
        f"Zero-gap = 0/{len(gap_res):,} (0%)"
    )
    ax1.text(0.97, 0.95, annotation, transform=ax1.transAxes,
             ha="right", va="top", fontsize=10,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="gray", alpha=0.9))

    ax1.legend(loc="upper center", fontsize=9, framealpha=0.9)
    sns.despine(ax=ax1)
    fig.tight_layout()
    save_fig(fig, "emergence_gap_distribution")
    plt.close(fig)


def figure2_kf_required_vs_empirical(kf_df):
    """scatter of kf_min vs segment size, colored by Galpern/validation."""
    print("figure 2: kf required vs empirical range")

    fig, ax = plt.subplots(figsize=(8, 5.5))

    # empirical band
    kf_sorted = kf_df.sort_values("seg_size")
    grouped = kf_sorted.groupby("seg_size").agg(
        emp_low=("empirical_kf_range_log10_low", "min"),
        emp_high=("empirical_kf_range_log10_high", "max"),
    ).reset_index()

    ax.fill_between(grouped["seg_size"], grouped["emp_low"], grouped["emp_high"],
                    color="lightgray", alpha=0.5, label="Empirical $k_f$ range", zorder=1)
    ax.plot(grouped["seg_size"], grouped["emp_low"], color="gray", linewidth=0.8, zorder=2)
    ax.plot(grouped["seg_size"], grouped["emp_high"], color="gray", linewidth=0.8, zorder=2)

    # scatter: split by Galpern vs validation
    galp = kf_df[kf_df["protein_id"].apply(is_galpern)]
    val = kf_df[~kf_df["protein_id"].apply(is_galpern)]

    ax.scatter(galp["seg_size"], galp["kf_min_log10"],
               color=GALPERN_COLOR, s=30, alpha=0.7, edgecolors="darkblue", linewidth=0.3,
               label=f"Galpern (n={len(galp)})", zorder=3)
    ax.scatter(val["seg_size"], val["kf_min_log10"],
               color=VALIDATION_COLOR, s=30, alpha=0.7, edgecolors="darkorange", linewidth=0.3,
               label=f"Validation (n={len(val)})", zorder=3)

    ax.set_xlabel("Foldon size (residues)")
    ax.set_ylabel(r"$\log_{10}(k_f$ / s$^{-1})$")

    nw = kf_df["within_empirical_range"].sum()
    nt = len(kf_df)
    ax.text(0.97, 0.05, f"{nw}/{nt} ({100*nw/nt:.0f}%) within empirical range",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=11,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="gray", alpha=0.9))

    ax.legend(loc="upper right", framealpha=0.9, fontsize=9)
    sns.despine(ax=ax)
    fig.tight_layout()
    save_fig(fig, "kf_required_vs_empirical")
    plt.close(fig)


def figure3_gap_vs_separation(gaps_df):
    """scatter of emergence gap vs segment separation."""
    print("figure 3: gap vs segment separation")

    df = gaps_df.copy()
    df["seg_separation"] = (df["seg_j"] - df["seg_i"]).abs()
    df["source"] = df["protein_id"].apply(lambda p: "Galpern" if is_galpern(p) else "Validation")

    fig, ax = plt.subplots(figsize=(8, 5.5))

    # separate by source
    for source, color, alpha in [("Galpern", GALPERN_COLOR, 0.5), ("Validation", VALIDATION_COLOR, 0.4)]:
        sub = df[df["source"] == source]
        ax.scatter(sub["seg_separation"], sub["emergence_gap_residues"],
                   color=color, s=15, alpha=alpha, edgecolors="none", label=source)

    # trend line (all data)
    x = df["seg_separation"].values
    y = df["emergence_gap_residues"].values
    slope, intercept, r_val, p_val, std_err = stats.linregress(x, y)
    x_fit = np.linspace(x.min(), x.max(), 100)
    y_fit = slope * x_fit + intercept
    ax.plot(x_fit, y_fit, color="black", linewidth=1.5, linestyle="--", zorder=5)

    ax.set_xlabel("Segment index separation")
    ax.set_ylabel("Emergence gap (residues)")

    p_str = f"{p_val:.2e}" if p_val > 0 else "< 1e-300"
    annotation = f"r = {r_val:.3f}, p {p_str}\nn = {len(df):,}"
    ax.text(0.03, 0.95, annotation, transform=ax.transAxes,
            ha="left", va="top", fontsize=11,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="gray", alpha=0.9))

    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
    sns.despine(ax=ax)
    fig.tight_layout()
    save_fig(fig, "gap_vs_separation")
    plt.close(fig)


def figure4_folding_timelines(seg_df):
    """timelines for 6 example proteins: 3 Galpern + 3 validation."""
    print("figure 4: folding timelines (6 examples)")

    TUNNEL_LEN = 30

    # pick 6 example proteins
    # 3 Galpern: F01 (small), F38 (ubiquitin, mostly B), F02 (large)
    # 3 validation: barnase (well-studied), R_9AVF (2 segments, tiny), R_5FAI (15 segs, large)
    candidates = ["F01", "F38", "F02", "barnase", "R_9AVF", "R_5FAI"]
    # only use proteins that exist in seg_df
    example_proteins = [p for p in candidates if p in seg_df["protein_id"].values]

    fig, axes = plt.subplots(len(example_proteins), 1, figsize=(10, 10))
    fig.subplots_adjust(hspace=0.55)

    for idx, pid in enumerate(example_proteins):
        ax = axes[idx]
        sub = seg_df[seg_df["protein_id"] == pid].sort_values("seg_idx")
        pdb_id = sub["pdb_id"].iloc[0]
        n_res = sub["n_res"].iloc[0]
        source = "Galpern" if is_galpern(pid) else "Validation"

        bar_y = 0.5
        bar_height = 0.35

        for _, row in sub.iterrows():
            seg_start = row["seg_start"]
            seg_end = row["seg_end"]
            seg_type = row["seg_type"]
            seg_idx = row["seg_idx"]

            color = TYPE_A_COLOR if seg_type == "A" else TYPE_B_COLOR
            ax.barh(bar_y, seg_end - seg_start, left=seg_start, height=bar_height,
                    color=color, edgecolor="white", linewidth=0.8, zorder=2)

            mid_x = (seg_start + seg_end) / 2
            ax.text(mid_x, bar_y, str(seg_idx), ha="center", va="center",
                    fontsize=8, fontweight="bold", color="white", zorder=3)

            emerge_x = seg_end + TUNNEL_LEN
            if emerge_x <= n_res + TUNNEL_LEN + 10:
                ax.axvline(emerge_x, color=color, linestyle=":", linewidth=0.8,
                           alpha=0.6, ymin=0.0, ymax=1.0, zorder=1)
                ax.plot(emerge_x, bar_y + bar_height / 2 + 0.08, marker="v",
                        color=color, markersize=5, zorder=4)

        ax.set_xlim(-5, max(n_res + TUNNEL_LEN + 15, n_res * 1.15))
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xlabel("Residue position" if idx == len(example_proteins) - 1 else "")
        ax.set_title(f"{pid}  ({pdb_id.upper()}, {n_res} res, {len(sub)} foldons, {source})",
                     fontsize=11, loc="left", fontweight="bold")
        sns.despine(ax=ax, left=True)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=TYPE_A_COLOR, edgecolor="white", label="Type A (autonomous)"),
        Patch(facecolor=TYPE_B_COLOR, edgecolor="white", label="Type B (cooperative)"),
    ]
    axes[0].legend(handles=legend_elements, loc="upper right", fontsize=9, framealpha=0.9)

    save_fig(fig, "folding_timelines")
    plt.close(fig)


def figure5_galpern_vs_validation(gaps_df):
    """comparison violin/box plot: Galpern vs validation emergence gaps."""
    print("figure 5: Galpern vs validation comparison")

    df = gaps_df.copy()
    df["source"] = df["protein_id"].apply(lambda p: "Galpern\n(n=37)" if is_galpern(p) else "Validation\n(n=67)")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    # left panel: emergence gap distributions
    galp = df[df["source"].str.startswith("Galpern")]["emergence_gap_residues"]
    val = df[df["source"].str.startswith("Validation")]["emergence_gap_residues"]

    parts = ax1.violinplot([galp, val], positions=[0, 1], showmedians=True, showextrema=False)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor([GALPERN_COLOR, VALIDATION_COLOR][i])
        pc.set_alpha(0.6)
    parts["cmedians"].set_color("black")

    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(["Galpern\n(n=37)", "Validation\n(n=67)"])
    ax1.set_ylabel("Emergence gap (residues)")
    ax1.set_title("Emergence gap distribution", fontsize=11, fontweight="bold")

    # medians
    ax1.text(0, galp.median() + 5, f"med={galp.median():.0f}", ha="center", fontsize=9)
    ax1.text(1, val.median() + 5, f"med={val.median():.0f}", ha="center", fontsize=9)

    # Mann-Whitney U test
    u_stat, u_pval = stats.mannwhitneyu(galp, val, alternative="two-sided")
    ax1.text(0.5, 0.95, f"Mann-Whitney p = {u_pval:.3f}", transform=ax1.transAxes,
             ha="center", va="top", fontsize=10,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", edgecolor="gray"))

    sns.despine(ax=ax1)

    # right panel: per-protein median gap
    per_protein = df.groupby("protein_id").agg(
        median_gap=("emergence_gap_residues", "median"),
        n_interfaces=("emergence_gap_residues", "count"),
    ).reset_index()
    per_protein["source"] = per_protein["protein_id"].apply(
        lambda p: "Galpern" if is_galpern(p) else "Validation"
    )

    galp_med = per_protein[per_protein["source"] == "Galpern"]["median_gap"]
    val_med = per_protein[per_protein["source"] == "Validation"]["median_gap"]

    ax2.hist(galp_med, bins=15, color=GALPERN_COLOR, alpha=0.6,
             label=f"Galpern (n={len(galp_med)})", edgecolor="white")
    ax2.hist(val_med, bins=15, color=VALIDATION_COLOR, alpha=0.6,
             label=f"Validation (n={len(val_med)})", edgecolor="white")

    ax2.set_xlabel("Per-protein median gap (residues)")
    ax2.set_ylabel("Number of proteins")
    ax2.set_title("Per-protein median emergence gap", fontsize=11, fontweight="bold")
    ax2.legend(fontsize=9, framealpha=0.9)

    u2, p2 = stats.mannwhitneyu(galp_med, val_med, alternative="two-sided")
    ax2.text(0.5, 0.95, f"Mann-Whitney p = {p2:.3f}", transform=ax2.transAxes,
             ha="center", va="top", fontsize=10,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", edgecolor="gray"))

    sns.despine(ax=ax2)
    fig.tight_layout()
    save_fig(fig, "galpern_vs_validation")
    plt.close(fig)


def main():
    print("loading data...")
    gaps_df = pd.read_csv(GAPS_CSV)
    kf_df = pd.read_csv(KF_CSV)
    seg_df = pd.read_csv(SEGTYPES_CSV)

    print(f"  emergence gaps: {len(gaps_df)} L3 interfaces")
    print(f"  kf required: {len(kf_df)} type A foldons")
    print(f"  segments: {len(seg_df)} across {seg_df['protein_id'].nunique()} proteins")
    print()

    figure1_emergence_gap_distribution(gaps_df)
    figure2_kf_required_vs_empirical(kf_df)
    figure3_gap_vs_separation(gaps_df)
    figure4_folding_timelines(seg_df)
    figure5_galpern_vs_validation(gaps_df)

    print()
    print(f"all figures saved to {FIG_DIR}/")


if __name__ == "__main__":
    main()
