"""
shared utilities for cotrans-layer pipeline.
codon table, amino acid mappings, data loaders.
"""

import pandas as pd
from collections import defaultdict
from pathlib import Path


# standard genetic code: codon → amino acid (or '*' for stop)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# reverse mapping: amino acid → list of synonymous codons (excludes stop)
AA_CODONS = defaultdict(list)
for _codon, _aa in CODON_TABLE.items():
    if _aa != '*':
        AA_CODONS[_aa].append(_codon)
AA_CODONS = dict(AA_CODONS)  # freeze as regular dict


def load_tai(path):
    """
    load per-codon tAI values from TSV.
    returns dict: codon → tAI (float)
    """
    df = pd.read_csv(path, sep='\t')
    return dict(zip(df['codon'], df['w_tai']))


def load_mu(path):
    """
    load per-codon mistranslation rates from TSV.
    returns dict: codon → mu (float)
    """
    df = pd.read_csv(path, sep='\t')
    return dict(zip(df['codon'], df['mu']))


def load_segment_types(path):
    """
    load segment_types.csv.
    columns: protein_id, pdb_id, n_res, seg_idx, seg_start, seg_end,
             seg_size, n_L2, n_L3, n_total, L2_ratio, E_L2, E_L3,
             seg_type, position, n_segments
    """
    return pd.read_csv(path)


def load_contacts(path, family_id=None):
    """
    load contacts_awsem.csv, optionally filter by family_id.
    columns: i, j, seq_sep, E_direct, E_water, E_protein, E_total,
             layer, seg_i, seg_j, family_id, pdb_id, n_res
    """
    df = pd.read_csv(path)
    if family_id is not None:
        df = df[df['family_id'] == family_id].copy()
    return df
