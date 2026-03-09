import pandas as pd
import numpy as np
import pytest
import sys
sys.path.insert(0, '/storage/kiran-stuff/codon_project/cotrans-layer')
from src.contact_analysis import extract_foldon_contacts, compute_all_foldon_cos

def test_extract_foldon_contacts_f01():
    """F01 segment 2 should have 22 L2 contacts (from segment_types)"""
    contacts_df = pd.read_csv('data/upstream/contacts_awsem.csv')
    seg_df = pd.read_csv('data/upstream/segment_types.csv')

    f01_contacts = contacts_df[contacts_df['family_id'] == 'F01']
    f01_segs = seg_df[seg_df['protein_id'] == 'F01']

    seg2 = f01_segs[f01_segs['seg_idx'] == 2].iloc[0]
    l2_pairs = extract_foldon_contacts(f01_contacts, seg_idx=2)
    assert len(l2_pairs) == seg2['n_L2']  # should be 22

def test_compute_all_foldon_cos_f01():
    """F01 should have CO for type A segments only (1, 2, 3), None for type B (0, 4)"""
    contacts_df = pd.read_csv('data/upstream/contacts_awsem.csv')
    seg_df = pd.read_csv('data/upstream/segment_types.csv')

    result = compute_all_foldon_cos('F01', contacts_df, seg_df)
    assert len(result) == 5  # 5 segments

    # type B segments (0, 4) should have co=None, kf=None
    assert result[0]['co'] is None
    assert result[4]['co'] is None

    # type A segments (1, 2, 3) should have co > 0 and kf > 0
    for idx in [1, 2, 3]:
        assert result[idx]['co'] is not None
        assert result[idx]['co'] > 0
        assert result[idx]['kf'] > 0
