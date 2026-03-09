"""
contact analysis for co-translational folding pipeline.

extracts per-foldon L2 contacts from the awsem contact data
and computes contact order / folding rate for each segment.
"""

import pandas as pd
from src.kinetic_model import compute_contact_order, compute_kf


def extract_foldon_contacts(contacts_df, seg_idx):
    """
    extract intra-foldon (L2) contacts for a given segment index.

    filters contacts_df to rows where:
      - layer == 'L2'
      - seg_i == seg_idx
      - seg_j == seg_idx
    (both residues must be within the same segment for L2)

    args:
        contacts_df: dataframe with columns i, j, layer, seg_i, seg_j
                     (already filtered to one protein)
        seg_idx: integer segment index to extract contacts for

    returns:
        list of (i, j) tuples — raw residue indices (0-indexed, protein-relative)
    """
    # filter to L2 contacts where both residues are in the target segment
    mask = (
        (contacts_df['layer'] == 'L2') &
        (contacts_df['seg_i'] == seg_idx) &
        (contacts_df['seg_j'] == seg_idx)
    )
    filtered = contacts_df[mask]

    # return as list of (i, j) tuples
    return list(zip(filtered['i'].values, filtered['j'].values))


def compute_all_foldon_cos(protein_id, contacts_df, seg_df):
    """
    compute contact order and folding rate for every segment of a protein.

    for type A segments (those with L2 contacts), computes CO using
    the absolute sequence separation |i - j| of L2 contacts and the
    segment size as the foldon length.

    for type B segments (no L2 contacts), returns co=None, kf=None.

    args:
        protein_id: string identifier (matches protein_id in seg_df,
                    family_id in contacts_df)
        contacts_df: full contacts dataframe (will be filtered by family_id)
        seg_df: full segment_types dataframe (will be filtered by protein_id)

    returns:
        list of dicts ordered by seg_idx, each containing:
          seg_idx, seg_type, co, kf, n_l2, seg_start, seg_end, seg_size
    """
    # filter to this protein
    prot_segs = seg_df[seg_df['protein_id'] == protein_id].sort_values('seg_idx')
    prot_contacts = contacts_df[contacts_df['family_id'] == protein_id]

    results = []
    for _, seg in prot_segs.iterrows():
        idx = seg['seg_idx']
        seg_size = seg['seg_size']

        # extract L2 contacts for this segment
        l2_contacts = extract_foldon_contacts(prot_contacts, seg_idx=idx)
        n_l2 = len(l2_contacts)

        if n_l2 == 0:
            # type B segment: no intra-foldon contacts
            co = None
            kf = None
        else:
            # type A segment: compute contact order and folding rate
            # contact i,j values are 0-indexed relative to protein (not foldon)
            # compute_contact_order uses |i - j| directly as sequence separation
            co = compute_contact_order(l2_contacts, foldon_length=seg_size)
            kf = compute_kf(co)

        results.append({
            'seg_idx': idx,
            'seg_type': seg['seg_type'],
            'co': co,
            'kf': kf,
            'n_l2': n_l2,
            'seg_start': seg['seg_start'],
            'seg_end': seg['seg_end'],
            'seg_size': seg_size,
        })

    return results
