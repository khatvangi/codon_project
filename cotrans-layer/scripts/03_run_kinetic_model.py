#!/usr/bin/env python3
"""
03_run_kinetic_model.py — run the O'Brien kinetic model for all proteins
under multiple parameter combinations (3 rates x 3 tunnels x 4 ku/kf ratios).

outputs:
  results/folding_curves/{protein_id}_cotrans.csv  (default params only)
  results/sensitivity/{protein_id}_sensitivity.csv  (all 36 combos)

usage:
  python scripts/03_run_kinetic_model.py --proteins F01,F38 --verbose
  python scripts/03_run_kinetic_model.py --all --workers 4
"""

import argparse
import sys
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from multiprocessing import Pool

BASE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE))

from src.kinetic_model import cotrans_folding_probability
from src.contact_analysis import compute_all_foldon_cos
from src.rate_computation import uniform_rates, tai_rates, shuffled_synonymous_rates
from src.cds_mapping import extract_codons
from src.utils import load_tai, load_contacts, load_segment_types

# ── paths ────────────────────────────────────────────────────────────
MANIFEST = BASE / "data" / "protein_manifest.csv"
CDS_DIR = BASE / "data" / "cds"
CONTACTS_PATH = BASE / "data" / "upstream" / "contacts_awsem.csv"
SEGTYPES_PATH = BASE / "data" / "upstream" / "segment_types.csv"
TAI_PATH = Path("/storage/kiran-stuff/codon_project/ecoli_tai_ws.tsv")

CURVES_DIR = BASE / "results" / "folding_curves"
SENS_DIR = BASE / "results" / "sensitivity"
CURVES_DIR.mkdir(parents=True, exist_ok=True)
SENS_DIR.mkdir(parents=True, exist_ok=True)

# ── parameter grid ───────────────────────────────────────────────────
RATE_REGIMES = ['uniform', 'tai', 'shuffled']
TUNNEL_LENGTHS = [25, 30, 40]
KU_KF_RATIOS = [0.0001, 0.001, 0.01, 0.1]
N_SHUFFLES = 100

# default params for folding curves output
DEFAULT_RATE = 'uniform'
DEFAULT_TUNNEL = 30
DEFAULT_KU_KF = 0.001


def read_fasta(fasta_path):
    """
    read a simple single-sequence FASTA file.
    returns the sequence string (no header, no newlines).
    """
    lines = fasta_path.read_text().strip().split('\n')
    # skip header line(s) starting with '>'
    seq_lines = [l.strip() for l in lines if not l.startswith('>')]
    return ''.join(seq_lines)


def compute_tau_array(codons, regime, tai_dict, seed=None):
    """
    compute per-codon elongation times for the given rate regime.
    returns numpy array of tau values.
    """
    if regime == 'uniform':
        return uniform_rates(codons)
    elif regime == 'tai':
        return tai_rates(codons, tai_dict)
    elif regime == 'shuffled':
        return shuffled_synonymous_rates(codons, tai_dict, seed=seed)
    else:
        raise ValueError(f"unknown rate regime: {regime}")


def run_one_protein(args):
    """
    process a single protein: compute folding curves and sensitivity table.
    args is a tuple: (protein_id, contacts_df, seg_df, tai_dict, verbose)
    """
    protein_id, contacts_df, seg_df, tai_dict, verbose = args

    # load CDS and extract codons
    fasta_path = CDS_DIR / f"{protein_id}.fasta"
    if not fasta_path.exists():
        print(f"  [{protein_id}] SKIP: no FASTA file")
        return None

    cds_seq = read_fasta(fasta_path)
    codons = extract_codons(cds_seq)
    n_codons = len(codons)

    if verbose:
        print(f"  [{protein_id}] {n_codons} codons")

    # compute per-foldon contact order and folding rates
    foldon_data = compute_all_foldon_cos(protein_id, contacts_df, seg_df)

    if not foldon_data:
        print(f"  [{protein_id}] SKIP: no segment data")
        return None

    if verbose:
        for fd in foldon_data:
            kf_str = f"{fd['kf']:.2f}" if fd['kf'] is not None else "None"
            print(f"    seg {fd['seg_idx']} ({fd['seg_type']}): "
                  f"n_l2={fd['n_l2']}, co={fd['co']}, kf={kf_str}")

    # ── folding curves (default parameters only) ─────────────────────
    tau_default = compute_tau_array(codons, DEFAULT_RATE, tai_dict)
    curves = {'chain_length': list(range(n_codons))}

    for fd in foldon_data:
        col = f"P_folded_seg_{fd['seg_idx']}"
        if fd['kf'] is None:
            # type B: no folding
            curves[col] = np.zeros(n_codons).tolist()
        else:
            ku = DEFAULT_KU_KF * fd['kf']
            P = cotrans_folding_probability(
                fd['seg_start'], fd['seg_end'], fd['kf'], ku,
                tau_default, tunnel_length=DEFAULT_TUNNEL
            )
            curves[col] = P.tolist()

    curves_df = pd.DataFrame(curves)
    curves_path = CURVES_DIR / f"{protein_id}_cotrans.csv"
    curves_df.to_csv(curves_path, index=False)

    if verbose:
        print(f"    saved folding curves: {curves_path.name}")

    # ── sensitivity analysis (all parameter combinations) ────────────
    sens_rows = []

    for regime in RATE_REGIMES:
        for tunnel in TUNNEL_LENGTHS:
            for ku_kf in KU_KF_RATIOS:
                # compute tau arrays (for shuffled: need multiple)
                if regime == 'shuffled':
                    # precompute all shuffled tau arrays
                    tau_shuffles = [
                        compute_tau_array(codons, 'shuffled', tai_dict, seed=s)
                        for s in range(N_SHUFFLES)
                    ]
                else:
                    tau_single = compute_tau_array(codons, regime, tai_dict)

                for fd in foldon_data:
                    seg_idx = fd['seg_idx']
                    seg_type = fd['seg_type']

                    if fd['kf'] is None:
                        # type B foldon: no folding possible
                        row = {
                            'protein_id': protein_id,
                            'seg_idx': seg_idx,
                            'seg_type': seg_type,
                            'n_l2': fd['n_l2'],
                            'co': fd['co'],
                            'kf': fd['kf'],
                            'rate_regime': regime,
                            'tunnel_length': tunnel,
                            'ku_kf_ratio': ku_kf,
                            'p_folded_at_emergence': 0.0,
                            'p_folded_final': 0.0,
                        }
                        if regime == 'shuffled':
                            row['p_folded_at_emergence_std'] = 0.0
                        sens_rows.append(row)
                        continue

                    # type A foldon: compute folding probability
                    ku = ku_kf * fd['kf']
                    emergence_idx = fd['seg_end'] + tunnel

                    if regime == 'shuffled':
                        # run N_SHUFFLES and average
                        p_emerge_values = []
                        p_final_values = []

                        for tau_shuf in tau_shuffles:
                            P = cotrans_folding_probability(
                                fd['seg_start'], fd['seg_end'],
                                fd['kf'], ku, tau_shuf,
                                tunnel_length=tunnel
                            )
                            # p_folded at emergence
                            if emergence_idx < n_codons:
                                p_emerge_values.append(P[emergence_idx])
                            else:
                                p_emerge_values.append(np.nan)
                            p_final_values.append(P[-1])

                        # suppress warnings for all-NaN slices (foldons that never emerge)
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", RuntimeWarning)
                            p_emerge_mean = np.nanmean(p_emerge_values)
                            p_emerge_std = np.nanstd(p_emerge_values)
                            p_final_mean = np.nanmean(p_final_values)

                        row = {
                            'protein_id': protein_id,
                            'seg_idx': seg_idx,
                            'seg_type': seg_type,
                            'n_l2': fd['n_l2'],
                            'co': co if (co := fd['co']) is not None else np.nan,
                            'kf': fd['kf'],
                            'rate_regime': regime,
                            'tunnel_length': tunnel,
                            'ku_kf_ratio': ku_kf,
                            'p_folded_at_emergence': p_emerge_mean,
                            'p_folded_final': p_final_mean,
                            'p_folded_at_emergence_std': p_emerge_std,
                        }
                    else:
                        # uniform or tai: single run
                        P = cotrans_folding_probability(
                            fd['seg_start'], fd['seg_end'],
                            fd['kf'], ku, tau_single,
                            tunnel_length=tunnel
                        )

                        if emergence_idx < n_codons:
                            p_emerge = P[emergence_idx]
                        else:
                            p_emerge = np.nan

                        row = {
                            'protein_id': protein_id,
                            'seg_idx': seg_idx,
                            'seg_type': seg_type,
                            'n_l2': fd['n_l2'],
                            'co': fd['co'] if fd['co'] is not None else np.nan,
                            'kf': fd['kf'],
                            'rate_regime': regime,
                            'tunnel_length': tunnel,
                            'ku_kf_ratio': ku_kf,
                            'p_folded_at_emergence': p_emerge,
                            'p_folded_final': P[-1],
                        }

                    sens_rows.append(row)

    sens_df = pd.DataFrame(sens_rows)
    sens_path = SENS_DIR / f"{protein_id}_sensitivity.csv"
    sens_df.to_csv(sens_path, index=False)

    if verbose:
        n_a = sum(1 for fd in foldon_data if fd['kf'] is not None)
        n_b = sum(1 for fd in foldon_data if fd['kf'] is None)
        print(f"    saved sensitivity: {sens_path.name} "
              f"({len(sens_rows)} rows, {n_a} type-A, {n_b} type-B foldons)")

    return protein_id


def main():
    parser = argparse.ArgumentParser(
        description="run O'Brien kinetic model for all proteins")
    parser.add_argument('--proteins', type=str, default=None,
                        help='comma-separated protein IDs (e.g. F01,F38)')
    parser.add_argument('--all', action='store_true',
                        help='run all proteins with cds_status=success')
    parser.add_argument('--verbose', action='store_true',
                        help='print per-foldon results')
    parser.add_argument('--workers', type=int, default=1,
                        help='number of parallel workers (default 1)')
    args = parser.parse_args()

    if not args.proteins and not args.all:
        parser.error("specify --proteins or --all")

    # load shared data
    manifest = pd.read_csv(MANIFEST)
    contacts_df = load_contacts(CONTACTS_PATH)
    seg_df = load_segment_types(SEGTYPES_PATH)
    tai_dict = load_tai(TAI_PATH)

    print(f"loaded: {len(contacts_df)} contacts, "
          f"{len(seg_df)} segment rows, "
          f"{len(tai_dict)} tAI values")

    # determine which proteins to run
    if args.all:
        protein_ids = manifest[manifest['cds_status'] == 'success']['protein_id'].tolist()
    else:
        protein_ids = [p.strip() for p in args.proteins.split(',')]

    print(f"running {len(protein_ids)} proteins: {', '.join(protein_ids)}\n")

    # build argument tuples
    job_args = [
        (pid, contacts_df, seg_df, tai_dict, args.verbose)
        for pid in protein_ids
    ]

    # run (serial or parallel)
    if args.workers > 1:
        with Pool(args.workers) as pool:
            results = pool.map(run_one_protein, job_args)
    else:
        results = [run_one_protein(a) for a in job_args]

    # summary
    completed = [r for r in results if r is not None]
    skipped = len(results) - len(completed)
    print(f"\n{'='*60}")
    print(f"done: {len(completed)} completed, {skipped} skipped")
    print(f"folding curves: {CURVES_DIR}")
    print(f"sensitivity:    {SENS_DIR}")


if __name__ == '__main__':
    main()
