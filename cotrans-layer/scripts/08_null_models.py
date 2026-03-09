#!/usr/bin/env python3
"""
08_null_models.py — null models for emergence gap analysis

tests whether real foldon segmentations produce systematically larger
emergence gaps than random segmentations. two null models:

  null A: random contiguous segmentation (same k, random cut points)
  null B: shuffled segment sizes (same sizes, random order)

for each null iteration, contacts are reclassified as intra-segment (L2)
or inter-segment (L3) based on the null segmentation, and emergence gaps
are recomputed from scratch.

output:
  results/summary_104/null_model_results.csv
  results/summary_104/null_model_distributions.csv
"""

import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations

# ── paths ────────────────────────────────────────────────────────────
BASE = Path("/storage/kiran-stuff/codon_project/cotrans-layer")
SEG_PATH = BASE / "data" / "extended_segment_types.csv"
CONTACTS_PATH = Path("/storage/kiran-stuff/foldon_project/layer_results/energy_decomposition.csv")
OUT_DIR = BASE / "results" / "summary_104"

N_ITER = 1000
SEED = 42


# ── helper functions ─────────────────────────────────────────────────

def get_segment_boundaries(seg_df, protein_id):
    """
    return sorted array of (start, end) tuples for a protein's segments.
    both start and end are 0-indexed inclusive.
    """
    prot = seg_df[seg_df["protein_id"] == protein_id].sort_values("seg_idx")
    boundaries = list(zip(prot["seg_start"].values, prot["seg_end"].values))
    return boundaries


def assign_contacts_to_segments(contact_i, contact_j, seg_ends):
    """
    vectorized segment assignment using searchsorted.

    seg_ends: sorted array of inclusive segment end positions.
    for a residue r, its segment index is the first seg_end >= r,
    i.e. np.searchsorted(seg_ends, r, side='left').

    returns (seg_i_arr, seg_j_arr) — segment indices for each contact.
    """
    seg_i = np.searchsorted(seg_ends, contact_i, side="left")
    seg_j = np.searchsorted(seg_ends, contact_j, side="left")
    return seg_i, seg_j


def compute_gap_metrics(seg_i_arr, seg_j_arr, seg_ends):
    """
    from segment assignments for all contacts, find unique L3 interface
    pairs and compute emergence gap metrics.

    returns dict with min_gap, median_gap, n_l3_interfaces, n_type_a.
    returns None if no L3 interfaces exist.
    """
    # find inter-segment contacts (L3 analogs)
    mask = seg_i_arr != seg_j_arr
    if not np.any(mask):
        return None

    # unique L3 interface pairs: (min_seg, max_seg)
    l3_i = seg_i_arr[mask]
    l3_j = seg_j_arr[mask]
    # canonicalize so smaller index first
    lo = np.minimum(l3_i, l3_j)
    hi = np.maximum(l3_i, l3_j)
    # unique pairs via set of tuples
    pairs = set(zip(lo.tolist(), hi.tolist()))

    if len(pairs) == 0:
        return None

    # compute emergence gaps: |seg_end[a] - seg_end[b]| for each pair
    gaps = np.array([abs(int(seg_ends[a]) - int(seg_ends[b])) for a, b in pairs])

    # count type A segments: segments with at least 1 intra-segment contact
    intra_mask = seg_i_arr == seg_j_arr
    if np.any(intra_mask):
        intra_segs = np.unique(seg_i_arr[intra_mask])
        n_type_a = len(intra_segs)
    else:
        n_type_a = 0

    return {
        "min_gap": int(gaps.min()),
        "median_gap": float(np.median(gaps)),
        "n_l3_interfaces": len(pairs),
        "n_type_a": n_type_a,
    }


def generate_null_a(n_res, k, rng):
    """
    null A: random contiguous segmentation.
    choose k-1 unique cut points from {1..n_res-1}, sort them.
    returns seg_ends array (0-indexed inclusive end positions).
    """
    if k == 1:
        return np.array([n_res - 1])
    # choose k-1 cut points from 1..n_res-1
    cuts = rng.choice(np.arange(1, n_res), size=k - 1, replace=False)
    cuts = np.sort(cuts)
    # segment ends: cut[0]-1, cut[1]-1, ..., n_res-1
    # segments are [0, cut[0]-1], [cut[0], cut[1]-1], ..., [cut[-1], n_res-1]
    # so seg_ends = [cut[0]-1, cut[1]-1, ..., n_res-1]
    seg_ends = np.empty(k, dtype=np.int64)
    seg_ends[:-1] = cuts - 1
    seg_ends[-1] = n_res - 1
    return seg_ends


def generate_null_b(sizes, rng):
    """
    null B: shuffled segment sizes.
    same multiset of sizes, randomly permuted order.
    returns seg_ends array (0-indexed inclusive end positions).
    """
    perm = rng.permutation(sizes)
    seg_ends = np.cumsum(perm) - 1
    return seg_ends


# ── main ─────────────────────────────────────────────────────────────

def main():
    print("loading data...")
    seg_df = pd.read_csv(SEG_PATH)
    contacts_df = pd.read_csv(CONTACTS_PATH)

    protein_ids = seg_df["protein_id"].unique()
    print(f"  {len(protein_ids)} proteins, {len(contacts_df)} contacts")

    rng = np.random.default_rng(SEED)

    # pre-index contacts per protein as numpy arrays
    contacts_by_protein = {}
    for pid, grp in contacts_df.groupby("protein_id"):
        contacts_by_protein[pid] = (
            grp["contact_i"].values.astype(np.int64),
            grp["contact_j"].values.astype(np.int64),
        )

    results_rows = []
    dist_rows = []

    print(f"running null models ({N_ITER} iterations each)...")
    print()

    for idx, pid in enumerate(protein_ids):
        if idx % 10 == 0:
            print(f"  protein {idx+1}/{len(protein_ids)}: {pid}")

        # get real segmentation
        boundaries = get_segment_boundaries(seg_df, pid)
        k = len(boundaries)
        n_res = seg_df[seg_df["protein_id"] == pid]["n_res"].iloc[0]

        # skip single-segment proteins (no L3 possible)
        if k <= 1:
            continue

        # real segment sizes
        sizes = np.array([end - start + 1 for start, end in boundaries])
        # real seg_ends
        real_seg_ends = np.array([end for _, end in boundaries])

        # get contacts for this protein
        if pid not in contacts_by_protein:
            continue
        ci, cj = contacts_by_protein[pid]
        if len(ci) == 0:
            continue

        # compute real metrics
        real_seg_i, real_seg_j = assign_contacts_to_segments(ci, cj, real_seg_ends)
        real_metrics = compute_gap_metrics(real_seg_i, real_seg_j, real_seg_ends)
        if real_metrics is None:
            continue

        # run null models
        null_a_min_gaps = np.empty(N_ITER)
        null_a_median_gaps = np.empty(N_ITER)
        null_a_n_l3 = np.empty(N_ITER, dtype=np.int64)
        null_a_n_type_a = np.empty(N_ITER, dtype=np.int64)

        null_b_min_gaps = np.empty(N_ITER)
        null_b_median_gaps = np.empty(N_ITER)
        null_b_n_l3 = np.empty(N_ITER, dtype=np.int64)
        null_b_n_type_a = np.empty(N_ITER, dtype=np.int64)

        for it in range(N_ITER):
            # ── null A ──
            se_a = generate_null_a(n_res, k, rng)
            si_a, sj_a = assign_contacts_to_segments(ci, cj, se_a)
            m_a = compute_gap_metrics(si_a, sj_a, se_a)
            if m_a is not None:
                null_a_min_gaps[it] = m_a["min_gap"]
                null_a_median_gaps[it] = m_a["median_gap"]
                null_a_n_l3[it] = m_a["n_l3_interfaces"]
                null_a_n_type_a[it] = m_a["n_type_a"]
            else:
                null_a_min_gaps[it] = np.nan
                null_a_median_gaps[it] = np.nan
                null_a_n_l3[it] = 0
                null_a_n_type_a[it] = 0

            # ── null B ──
            se_b = generate_null_b(sizes, rng)
            si_b, sj_b = assign_contacts_to_segments(ci, cj, se_b)
            m_b = compute_gap_metrics(si_b, sj_b, se_b)
            if m_b is not None:
                null_b_min_gaps[it] = m_b["min_gap"]
                null_b_median_gaps[it] = m_b["median_gap"]
                null_b_n_l3[it] = m_b["n_l3_interfaces"]
                null_b_n_type_a[it] = m_b["n_type_a"]
            else:
                null_b_min_gaps[it] = np.nan
                null_b_median_gaps[it] = np.nan
                null_b_n_l3[it] = 0
                null_b_n_type_a[it] = 0

        # p-values: fraction of null iterations with min_gap >= real min_gap
        # (one-sided test: is the real gap unusually large?)
        valid_a = ~np.isnan(null_a_min_gaps)
        valid_b = ~np.isnan(null_b_min_gaps)

        if valid_a.sum() > 0:
            pval_a_min = (null_a_min_gaps[valid_a] >= real_metrics["min_gap"]).mean()
            pval_a_median = (null_a_median_gaps[valid_a] >= real_metrics["median_gap"]).mean()
        else:
            pval_a_min = np.nan
            pval_a_median = np.nan

        if valid_b.sum() > 0:
            pval_b_min = (null_b_min_gaps[valid_b] >= real_metrics["min_gap"]).mean()
            pval_b_median = (null_b_median_gaps[valid_b] >= real_metrics["median_gap"]).mean()
        else:
            pval_b_min = np.nan
            pval_b_median = np.nan

        results_rows.append({
            "protein_id": pid,
            "n_res": n_res,
            "n_segments": k,
            "real_min_gap": real_metrics["min_gap"],
            "real_median_gap": real_metrics["median_gap"],
            "real_n_l3": real_metrics["n_l3_interfaces"],
            "real_n_type_a": real_metrics["n_type_a"],
            "null_a_pvalue_min_gap": pval_a_min,
            "null_a_pvalue_median_gap": pval_a_median,
            "null_a_median_min_gap": np.nanmedian(null_a_min_gaps),
            "null_a_mean_min_gap": np.nanmean(null_a_min_gaps),
            "null_b_pvalue_min_gap": pval_b_min,
            "null_b_pvalue_median_gap": pval_b_median,
            "null_b_median_min_gap": np.nanmedian(null_b_min_gaps),
            "null_b_mean_min_gap": np.nanmean(null_b_min_gaps),
        })

        # save summary percentiles for distribution file
        for label, mg, medg in [
            ("A", null_a_min_gaps, null_a_median_gaps),
            ("B", null_b_min_gaps, null_b_median_gaps),
        ]:
            valid = ~np.isnan(mg)
            if valid.sum() == 0:
                continue
            dist_rows.append({
                "protein_id": pid,
                "null_type": label,
                "min_gap_p05": np.nanpercentile(mg, 5),
                "min_gap_p25": np.nanpercentile(mg, 25),
                "min_gap_p50": np.nanpercentile(mg, 50),
                "min_gap_p75": np.nanpercentile(mg, 75),
                "min_gap_p95": np.nanpercentile(mg, 95),
                "median_gap_p50": np.nanpercentile(medg, 50),
            })

    # ── save results ─────────────────────────────────────────────────
    results_df = pd.DataFrame(results_rows)
    dist_df = pd.DataFrame(dist_rows)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    results_path = OUT_DIR / "null_model_results.csv"
    results_df.to_csv(results_path, index=False)
    print(f"\nsaved: {results_path}")

    dist_path = OUT_DIR / "null_model_distributions.csv"
    dist_df.to_csv(dist_path, index=False)
    print(f"saved: {dist_path}")

    # ── summary ──────────────────────────────────────────────────────
    print()
    print("=" * 70)
    print("NULL MODEL RESULTS SUMMARY")
    print("=" * 70)
    print()

    n_proteins = len(results_df)
    print(f"proteins analyzed: {n_proteins}")
    print()

    # ── null A summary ───────────────────────────────────────────────
    print("── null A: random contiguous segmentation ──")
    valid_a = results_df.dropna(subset=["null_a_pvalue_min_gap"])
    n_real_gt_null_a = (valid_a["real_min_gap"] > valid_a["null_a_median_min_gap"]).sum()
    print(f"  proteins with real min_gap > null median: {n_real_gt_null_a}/{len(valid_a)} "
          f"({100*n_real_gt_null_a/len(valid_a):.1f}%)")

    sig_a = (valid_a["null_a_pvalue_min_gap"] < 0.05).sum()
    print(f"  proteins with p < 0.05 (min_gap): {sig_a}/{len(valid_a)} "
          f"({100*sig_a/len(valid_a):.1f}%)")

    sig_a_median = (valid_a["null_a_pvalue_median_gap"] < 0.05).sum()
    print(f"  proteins with p < 0.05 (median_gap): {sig_a_median}/{len(valid_a)} "
          f"({100*sig_a_median/len(valid_a):.1f}%)")

    # fisher's combined p-value for null A (min_gap)
    from scipy.stats import combine_pvalues
    pvals_a = valid_a["null_a_pvalue_min_gap"].values
    # clamp to avoid log(0)
    pvals_a = np.clip(pvals_a, 1e-10, 1.0)
    fisher_stat_a, fisher_p_a = combine_pvalues(pvals_a, method="fisher")
    print(f"  fisher combined p-value (min_gap): {fisher_p_a:.2e}")
    print()

    # ── null B summary ───────────────────────────────────────────────
    print("── null B: shuffled segment sizes ──")
    valid_b = results_df.dropna(subset=["null_b_pvalue_min_gap"])
    n_real_gt_null_b = (valid_b["real_min_gap"] > valid_b["null_b_median_min_gap"]).sum()
    print(f"  proteins with real min_gap > null median: {n_real_gt_null_b}/{len(valid_b)} "
          f"({100*n_real_gt_null_b/len(valid_b):.1f}%)")

    sig_b = (valid_b["null_b_pvalue_min_gap"] < 0.05).sum()
    print(f"  proteins with p < 0.05 (min_gap): {sig_b}/{len(valid_b)} "
          f"({100*sig_b/len(valid_b):.1f}%)")

    sig_b_median = (valid_b["null_b_pvalue_median_gap"] < 0.05).sum()
    print(f"  proteins with p < 0.05 (median_gap): {sig_b_median}/{len(valid_b)} "
          f"({100*sig_b_median/len(valid_b):.1f}%)")

    pvals_b = valid_b["null_b_pvalue_min_gap"].values
    pvals_b = np.clip(pvals_b, 1e-10, 1.0)
    fisher_stat_b, fisher_p_b = combine_pvalues(pvals_b, method="fisher")
    print(f"  fisher combined p-value (min_gap): {fisher_p_b:.2e}")
    print()

    # ── aggregate comparison ─────────────────────────────────────────
    print("── aggregate ──")
    print(f"  real min_gap:   mean={results_df['real_min_gap'].mean():.1f}, "
          f"median={results_df['real_min_gap'].median():.1f}")
    print(f"  null A min_gap: mean={results_df['null_a_mean_min_gap'].mean():.1f}, "
          f"median={results_df['null_a_median_min_gap'].median():.1f}")
    print(f"  null B min_gap: mean={results_df['null_b_mean_min_gap'].mean():.1f}, "
          f"median={results_df['null_b_median_min_gap'].median():.1f}")
    print()

    # ── per-protein p-value table (first 20) ─────────────────────────
    print("── per-protein p-values (first 20) ──")
    print(f"{'protein':<20} {'n_seg':>5} {'real_min':>8} {'pval_A':>8} {'pval_B':>8}")
    print("-" * 55)
    for _, row in results_df.head(20).iterrows():
        print(f"{row['protein_id']:<20} {row['n_segments']:>5.0f} "
              f"{row['real_min_gap']:>8.0f} "
              f"{row['null_a_pvalue_min_gap']:>8.3f} "
              f"{row['null_b_pvalue_min_gap']:>8.3f}")
    if len(results_df) > 20:
        print(f"  ... ({len(results_df) - 20} more proteins)")
    print()

    # ── the key question ─────────────────────────────────────────────
    print("=" * 70)
    print("KEY QUESTION: do real segmentations produce larger minimum gaps")
    print("than random segmentations?")
    print()
    if fisher_p_a < 0.05 and n_real_gt_null_a > len(valid_a) / 2:
        print("  null A: YES — real segmentations have significantly larger gaps")
        print(f"          than random k-segment partitions (fisher p = {fisher_p_a:.2e})")
    else:
        print("  null A: NO — real segmentations are not significantly different")
        print(f"          from random k-segment partitions (fisher p = {fisher_p_a:.2e})")

    if fisher_p_b < 0.05 and n_real_gt_null_b > len(valid_b) / 2:
        print("  null B: YES — the specific ordering of segment sizes matters")
        print(f"          (fisher p = {fisher_p_b:.2e})")
    else:
        print("  null B: NO — shuffling segment sizes does not significantly")
        print(f"          reduce gaps (fisher p = {fisher_p_b:.2e})")
    print("=" * 70)


if __name__ == "__main__":
    main()
