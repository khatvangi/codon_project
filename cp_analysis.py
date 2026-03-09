#!/usr/bin/env python3
"""
cp_analysis.py — circular permutation layer analysis for protein S6 (1RIS).

tests the central prediction of the foldon-geometry framework:
the DP-optimal modular decomposition is an emergent consequence of
which residues are contiguous in the linear chain. circular permutation
changes contiguity without changing the 3D structure, so the same
contact energies get partitioned into different modules.

approach:
  1. load WT AWsemData (precomputed AWSEM weights from native structure)
  2. for each CP variant, create CP AWsemData by permuting coords/sequence
     and recomputing pair_mask (which encodes polymer connectivity)
  3. run standard delta_f scan + DP on CP data (unmodified pipeline)
  4. classify contacts using CP segments, compare with WT classification

the permutation trick: permuting the AWsemData arrays to CP coordinates
means the standard scan loop works unchanged. theta/sigma are distance-
dependent (invariant under permutation), pair_mask is recomputed for CP
connectivity, and burial weights are recomputed with CP |i-j| >= 2 filter.
"""

import sys
import os
import numpy as np
import pandas as pd
import pickle
from pathlib import Path
from collections import defaultdict

# add foldon project to path for LayerScan imports
FOLDON_DIR = Path("/storage/kiran-stuff/foldon_project")
sys.path.insert(0, str(FOLDON_DIR))
sys.path.insert(0, str(FOLDON_DIR / "layerscan"))

from layerscan.models import AWsemData, ScanOptions
from layerscan.scan import scan_delta_f_landscape, find_foldons
from layerscan.classify import compute_contacts, classify
from awsem_delta_f_scan import precompute_contact_weights, AA_MAP


# ============================================================
# CP-SPECIFIC FUNCTIONS
# ============================================================

def make_cp_data(wt_data: AWsemData, cut_site: int) -> AWsemData:
    """
    create AWsemData for a circular permutant.

    the CP chain starts at WT residue `cut_site` and wraps around:
    CP sequence: WT[cut_site], WT[cut_site+1], ..., WT[N-1], WT[0], ..., WT[cut_site-1]

    all 3D-structure-dependent quantities (theta, sigma) are permuted.
    pair_mask is recomputed for CP backbone connectivity.
    burial weights are recomputed with CP |i-j| >= 2 exclusion.

    args:
        wt_data: wild-type AWsemData
        cut_site: first residue of the CP chain (0-indexed, in WT numbering)

    returns:
        AWsemData with permuted arrays and fresh pair_mask
    """
    n = wt_data.n_residues
    # permutation: cp_pos p → wt_pos (p + cut_site) % n
    perm = np.array([(p + cut_site) % n for p in range(n)])

    # permute sequence and coords
    cp_seq = "".join(wt_data.sequence[i] for i in perm)
    cp_coords = wt_data.cb_coords[perm]
    cp_seq_types = wt_data.seq_types[perm]

    # recompute all distance-dependent weights from permuted coordinates.
    # this gives us:
    #   - theta_ij, thetaII_ij, sigma_w, sigma_p: same values, permuted indices
    #   - pair_mask: fresh, based on CP backbone connectivity (|i-j| >= 10)
    #   - burial_weights: recomputed with CP |i-j| >= 2 exclusion
    theta_ij, thetaII_ij, sigma_w, sigma_p, burial_w, pair_mask = \
        precompute_contact_weights(cp_coords, n)

    return AWsemData(
        sequence=cp_seq,
        n_residues=n,
        cb_coords=cp_coords,
        seq_types=cp_seq_types,
        gamma_d=wt_data.gamma_d,     # 20x20 matrices — protein-independent
        gamma_w=wt_data.gamma_w,
        gamma_p=wt_data.gamma_p,
        burial_gamma=wt_data.burial_gamma,
        theta_ij=theta_ij,
        thetaII_ij=thetaII_ij,
        sigma_w_ij=sigma_w,
        sigma_p_ij=sigma_p,
        burial_weights=burial_w,
        pair_mask=pair_mask,
    )


def cp_to_wt(cp_pos: int, cut_site: int, n: int) -> int:
    """convert CP position to WT residue index."""
    return (cp_pos + cut_site) % n


def wt_to_cp(wt_pos: int, cut_site: int, n: int) -> int:
    """convert WT residue index to CP position."""
    return (wt_pos - cut_site) % n


def map_segments_to_wt(cp_segments, cut_site, n):
    """
    map CP segment boundaries back to WT residue numbers.

    returns list of (cp_start, cp_end, wt_residues_set, wt_range_str)
    for human-readable comparison.
    """
    mapped = []
    for s, e in cp_segments:
        wt_residues = [cp_to_wt(p, cut_site, n) for p in range(s, e + 1)]
        # format as WT ranges (may wrap)
        wt_sorted = sorted(wt_residues)
        # detect if contiguous in WT
        ranges = []
        start = wt_sorted[0]
        prev = start
        for r in wt_sorted[1:]:
            if r == prev + 1:
                prev = r
            else:
                ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
                start = r
                prev = r
        ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
        wt_str = " + ".join(ranges)

        mapped.append({
            'cp_start': s, 'cp_end': e, 'cp_size': e - s + 1,
            'wt_residues': set(wt_residues),
            'wt_range': wt_str,
            'wraps': any(wt_residues[i] > wt_residues[i+1]
                         for i in range(len(wt_residues)-1)),
        })
    return mapped


def compare_layer_assignments(wt_contacts, cp_contacts, wt_segments,
                              cp_segments, cut_site, n):
    """
    compare L2/L3 assignments between WT and CP for the same contacts.

    contacts are identified by (wt_i, wt_j) pairs. in the CP contacts_df,
    the indices are in CP coordinates, so we convert back.

    returns:
        DataFrame with columns: wt_i, wt_j, E_direct, wt_layer, cp_layer, switched
    """
    # build WT lookup: (i, j) → layer, E_direct
    wt_lookup = {}
    for _, row in wt_contacts.iterrows():
        key = (int(row['i']), int(row['j']))
        wt_lookup[key] = {'layer': row['layer'], 'E_direct': row['E_direct'],
                          'E_total': row['E_total']}

    # build CP lookup: convert CP indices to WT, then look up
    cp_lookup = {}
    for _, row in cp_contacts.iterrows():
        wt_i = cp_to_wt(int(row['i']), cut_site, n)
        wt_j = cp_to_wt(int(row['j']), cut_site, n)
        # normalize order (i < j in WT numbering)
        if wt_i > wt_j:
            wt_i, wt_j = wt_j, wt_i
        cp_lookup[(wt_i, wt_j)] = {'layer': row['layer'], 'E_direct': row['E_direct']}

    # merge: contacts that exist in both WT and CP
    # note: some contacts may exist in one but not the other due to pair_mask change
    all_keys = set(wt_lookup.keys()) | set(cp_lookup.keys())

    rows = []
    for wt_i, wt_j in sorted(all_keys):
        wt_info = wt_lookup.get((wt_i, wt_j), {})
        cp_info = cp_lookup.get((wt_i, wt_j), {})

        wt_layer = wt_info.get('layer', 'absent')
        cp_layer = cp_info.get('layer', 'absent')
        e_direct = wt_info.get('E_direct', cp_info.get('E_direct', np.nan))
        e_total = wt_info.get('E_total', np.nan)

        rows.append({
            'wt_i': wt_i, 'wt_j': wt_j,
            'E_direct': e_direct, 'E_total': e_total,
            'wt_layer': wt_layer, 'cp_layer': cp_layer,
            'switched': wt_layer != cp_layer,
        })

    return pd.DataFrame(rows)


def run_cp_analysis(wt_data, wt_result, cut_site, label="CP"):
    """
    full CP analysis for one circular permutant.

    args:
        wt_data: AWsemData for wild-type
        wt_result: LayerScanResult for wild-type (for comparison)
        cut_site: CP cut position (0-indexed, WT numbering)
        label: human-readable label

    returns:
        dict with all CP results + comparison
    """
    n = wt_data.n_residues
    print(f"\n{'='*60}")
    print(f"  {label}: cut at WT residue {cut_site}")
    print(f"  CP chain: WT[{cut_site}..{n-1}] + WT[0..{cut_site-1}]")
    print(f"{'='*60}")

    # step 1: create CP AWsemData
    print("creating CP AWsemData...")
    cp_data = make_cp_data(wt_data, cut_site)
    print(f"  CP sequence: {cp_data.sequence[:20]}...{cp_data.sequence[-10:]}")
    print(f"  CP pair_mask contacts: {cp_data.pair_mask.sum() // 2}")

    # step 2: scan delta_f landscape
    print("scanning CP δf landscape...")
    options = ScanOptions()
    max_size = min(int(options.max_frag_frac * n), options.dp_max_size)
    delta_f, f_indep, f_context = scan_delta_f_landscape(cp_data, options)

    valid_entries = np.sum(~np.isnan(delta_f))
    print(f"  scanned: {valid_entries} valid δf entries")

    # step 3: DP-optimal segmentation
    print("finding CP-optimal foldons...")
    cp_segments, cp_cost = find_foldons(delta_f, n,
                                         min_size=options.dp_min_size,
                                         max_size=max_size)
    print(f"  {len(cp_segments)} foldons, total cost = {cp_cost:.4f}")

    # map segments to WT numbering for comparison
    mapped = map_segments_to_wt(cp_segments, cut_site, n)
    for i, m in enumerate(mapped):
        wraps_str = " [WRAPS]" if m['wraps'] else ""
        print(f"  CP-F{i}: CP[{m['cp_start']}-{m['cp_end']}] "
              f"= WT[{m['wt_range']}] size={m['cp_size']}{wraps_str}")

    # step 4: classify contacts
    print("classifying CP contacts...")
    cp_contacts = compute_contacts(cp_data)
    cp_contacts = classify(cp_contacts, cp_segments, n)

    n_l2 = len(cp_contacts[cp_contacts['layer'] == 'L2'])
    n_l3 = len(cp_contacts[cp_contacts['layer'] == 'L3'])
    print(f"  CP contacts: L2={n_l2}, L3={n_l3}")

    if n_l2 > 0:
        cp_l2_ed = cp_contacts[cp_contacts['layer'] == 'L2']['E_direct'].mean()
        cp_l3_ed = cp_contacts[cp_contacts['layer'] == 'L3']['E_direct'].mean()
        print(f"  CP L2 mean E_direct: {cp_l2_ed:.4f}")
        print(f"  CP L3 mean E_direct: {cp_l3_ed:.4f}")

    # step 5: compare with WT
    print("\ncomparing with WT...")
    wt_segments = [(f.start, f.end) for f in wt_result.foldons]
    comparison = compare_layer_assignments(
        wt_result.contacts_df, cp_contacts,
        wt_segments, cp_segments, cut_site, n)

    n_switched = comparison['switched'].sum()
    n_total = len(comparison)
    print(f"  total contacts in union: {n_total}")
    print(f"  contacts that switched layer: {n_switched} ({100*n_switched/n_total:.1f}%)")

    # break down switching by type
    for from_layer in ['L2', 'L3', 'gap', 'absent']:
        for to_layer in ['L2', 'L3', 'gap', 'absent']:
            if from_layer == to_layer:
                continue
            mask = (comparison['wt_layer'] == from_layer) & \
                   (comparison['cp_layer'] == to_layer)
            count = mask.sum()
            if count > 0:
                sub = comparison[mask]
                mean_ed = sub['E_direct'].mean()
                print(f"    {from_layer} → {to_layer}: {count} contacts, "
                      f"mean E_direct = {mean_ed:.4f}")

    # analyze contacts from WT-F2 (the module CP54 disrupts)
    # WT F2 = [44-71]
    f2_wt = set(range(44, 72))
    f2_contacts = comparison[
        comparison['wt_i'].isin(f2_wt) | comparison['wt_j'].isin(f2_wt)
    ]
    f2_switched = f2_contacts[f2_contacts['switched']]
    if len(f2_contacts) > 0:
        print(f"\n  WT-F2 [44-71] contacts: {len(f2_contacts)}")
        print(f"  WT-F2 switched: {len(f2_switched)} ({100*len(f2_switched)/len(f2_contacts):.1f}%)")
        if len(f2_switched) > 0:
            print(f"  WT-F2 switched mean E_direct: {f2_switched['E_direct'].mean():.4f}")

    return {
        'label': label,
        'cut_site': cut_site,
        'cp_data': cp_data,
        'cp_segments': cp_segments,
        'cp_cost': cp_cost,
        'cp_contacts': cp_contacts,
        'mapped_segments': mapped,
        'comparison': comparison,
        'delta_f': delta_f,
    }


# ============================================================
# MAIN
# ============================================================

def main():
    # paths
    wt_data_path = FOLDON_DIR / "cp_analysis" / "1ris_results" / "1ris.npz"
    wt_result_path = FOLDON_DIR / "cp_analysis" / "1ris_results" / "1ris" / "wt_result.pkl"
    output_dir = Path("/storage/kiran-stuff/codon_project/cp_results")
    output_dir.mkdir(parents=True, exist_ok=True)

    # load WT data
    print("loading WT AWsemData...")
    wt_data = AWsemData.load(str(wt_data_path))
    print(f"  {wt_data.n_residues} residues: {wt_data.sequence}")

    print("loading WT LayerScan result...")
    with open(wt_result_path, 'rb') as f:
        wt_result = pickle.load(f)

    print(f"\n=== WT segments ===")
    for fld in wt_result.foldons:
        print(f"  F{fld.index}: [{fld.start}-{fld.end}] size={fld.size} "
              f"type={fld.type_ab} n_L2={fld.n_l2} E_direct={fld.mean_e_direct}")

    # define CP variants to test
    # CP54: strongest experimental phenotype (cuts through WT-F2 [44-71])
    # CP13: cuts in WT-F1 [10-43] near the start
    # CP68: cuts near end of WT-F2 [44-71]
    cp_variants = [
        (13, "CP13 (within F1 nucleus)"),
        (68, "CP68 (near end of F2)"),
    ]

    results = {}
    for cut_site, label in cp_variants:
        res = run_cp_analysis(wt_data, wt_result, cut_site, label)
        results[cut_site] = res

        # save CP contacts and comparison
        res['cp_contacts'].to_csv(
            output_dir / f"cp{cut_site}_contacts.csv", index=False)
        res['comparison'].to_csv(
            output_dir / f"cp{cut_site}_comparison.csv", index=False)
        np.save(output_dir / f"cp{cut_site}_delta_f.npy", res['delta_f'])

    # save WT contacts for reference
    wt_result.contacts_df.to_csv(output_dir / "wt_contacts.csv", index=False)

    # print summary
    print(f"\n{'='*60}")
    print("  SUMMARY")
    print(f"{'='*60}")
    print(f"\nWT segments: {[(f.start, f.end) for f in wt_result.foldons]}")
    for cut_site, res in results.items():
        print(f"\n{res['label']}:")
        print(f"  segments (CP coords): {res['cp_segments']}")
        for m in res['mapped_segments']:
            print(f"    CP[{m['cp_start']}-{m['cp_end']}] = WT[{m['wt_range']}]")
        comp = res['comparison']
        n_sw = comp['switched'].sum()
        print(f"  switched: {n_sw}/{len(comp)} contacts ({100*n_sw/len(comp):.1f}%)")

    print(f"\nresults saved to: {output_dir}")


if __name__ == "__main__":
    main()
