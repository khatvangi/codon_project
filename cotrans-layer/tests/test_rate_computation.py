import numpy as np
import pytest
import sys
sys.path.insert(0, '/storage/kiran-stuff/codon_project/cotrans-layer')
from src.rate_computation import uniform_rates, tai_rates, shuffled_synonymous_rates

def test_uniform_rates():
    """all codons get tau=0.06s"""
    codons = ['ATG', 'AAA', 'GCT', 'TGA']
    rates = uniform_rates(codons, tau_mean=0.06)
    assert len(rates) == 4
    assert np.allclose(rates, 0.06)

def test_tai_rates_mean():
    """mean of tAI rates should equal tau_mean"""
    codons = ['ATG', 'AAA', 'GCT', 'GAT', 'TTT'] * 20
    tai_dict = {'ATG': 0.22, 'AAA': 1.0, 'GCT': 0.197, 'GAT': 0.295, 'TTT': 0.197}
    rates = tai_rates(codons, tai_dict, tau_mean=0.06, tai_floor=0.01)
    assert abs(np.mean(rates) - 0.06) < 0.001

def test_tai_rates_inverse_relationship():
    """higher tAI → lower tau (faster translation)"""
    codons = ['AAA', 'TTT']
    tai_dict = {'AAA': 1.0, 'TTT': 0.197}
    rates = tai_rates(codons, tai_dict, tau_mean=0.06, tai_floor=0.01)
    assert rates[0] < rates[1]  # AAA faster (lower tau)

def test_shuffled_preserves_aa():
    """shuffled synonymous should preserve amino acid identity"""
    codons = ['GCT', 'GCC', 'GCA', 'GCG']  # all Ala
    tai_dict = {'GCT': 0.197, 'GCC': 0.333, 'GCA': 0.5, 'GCG': 0.16}
    rates = shuffled_synonymous_rates(codons, tai_dict, tau_mean=0.06,
                                       tai_floor=0.01, seed=42)
    assert len(rates) == 4
    assert abs(np.mean(rates) - 0.06) < 0.01

def test_tai_floor_clips_extreme():
    """CGA has tAI ~0.00007 — floor should prevent extreme tau"""
    codons = ['CGA', 'AAA']
    tai_dict = {'CGA': 0.000067, 'AAA': 1.0}
    rates = tai_rates(codons, tai_dict, tau_mean=0.06, tai_floor=0.01)
    # without floor: CGA tau would be ~857s. with floor: manageable
    assert rates[0] < 10.0  # not absurdly slow
    assert rates[0] > rates[1]  # CGA still slower than AAA
