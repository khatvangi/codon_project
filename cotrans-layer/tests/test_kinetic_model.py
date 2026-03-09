import numpy as np
import pytest
import sys
sys.path.insert(0, '/storage/kiran-stuff/codon_project/cotrans-layer')
from src.kinetic_model import compute_contact_order, compute_kf, cotrans_folding_probability

def test_contact_order_simple():
    """two contacts at known separations in a 10-residue foldon"""
    contacts = [(0, 5), (2, 8)]  # seq_sep = 5, 6
    co = compute_contact_order(contacts, foldon_length=10)
    expected = (1.0 / (10 * 2)) * (5 + 6)  # = 0.55
    assert abs(co - expected) < 1e-10

def test_contact_order_empty():
    """zero contacts → CO undefined, return None"""
    co = compute_contact_order([], foldon_length=10)
    assert co is None

def test_kf_from_plaxco():
    """known CO → known kf via Plaxco relation"""
    co = 0.2
    kf = compute_kf(co)
    # ln(kf) = -17.0 * 0.2 + 8.0 = 4.6
    expected = np.exp(4.6)
    assert abs(kf - expected) < 1e-6

def test_pfold_before_emergence():
    """chain length before emergence → P_folded = 0"""
    codon_times = np.full(100, 0.06)
    P = cotrans_folding_probability(
        seg_start=10, seg_end=30, kf=100.0, ku=0.1,
        codon_times=codon_times, tunnel_length=30
    )
    # emergence at 30 + 30 = 60
    assert np.all(P[:61] == 0.0)

def test_pfold_approaches_equilibrium():
    """long time after emergence → P_folded → kf/(kf+ku)"""
    codon_times = np.full(200, 0.06)
    kf, ku = 100.0, 0.1
    P = cotrans_folding_probability(
        seg_start=10, seg_end=30, kf=kf, ku=ku,
        codon_times=codon_times, tunnel_length=30
    )
    P_eq = kf / (kf + ku)
    assert abs(P[-1] - P_eq) < 0.01

def test_pfold_monotonically_increasing():
    """P_folded should never decrease with chain length"""
    codon_times = np.full(150, 0.06)
    P = cotrans_folding_probability(
        seg_start=5, seg_end=20, kf=50.0, ku=0.05,
        codon_times=codon_times, tunnel_length=30
    )
    diffs = np.diff(P)
    assert np.all(diffs >= -1e-15)
