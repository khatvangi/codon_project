"""
rate computation for cotranslational folding simulations.
three models: uniform, tAI-based, and shuffled-synonymous control.
"""

import numpy as np
from src.utils import CODON_TABLE, AA_CODONS


def uniform_rates(codons, tau_mean=0.06):
    """
    assign identical elongation time to every codon.
    returns array of length len(codons), all equal to tau_mean.
    """
    return np.full(len(codons), tau_mean)


def tai_rates(codons, tai_dict, tau_mean=0.06, tai_floor=0.01):
    """
    compute per-codon elongation times from tAI values.

    logic:
      - look up tAI for each codon (missing codons get median tAI)
      - clip tAI to tai_floor to prevent extreme tau from near-zero tAI
        (e.g. CGA has tAI ~0.00007, which would give tau ~14000x the mean)
      - raw_tau = 1 / clipped_tAI  (inverse: high tAI = fast = low tau)
      - normalize so mean(tau) == tau_mean

    parameters:
      codons    : list of codon strings (e.g. ['ATG', 'AAA', ...])
      tai_dict  : dict codon -> tAI value
      tau_mean  : desired mean elongation time (seconds)
      tai_floor : minimum tAI value to prevent extreme outliers

    returns:
      numpy array of elongation times, one per codon
    """
    # fallback for codons not in tai_dict: use median of known values
    median_tai = np.median(list(tai_dict.values()))

    # look up tAI for each codon, clip to floor
    tai_values = np.array([
        max(tai_dict.get(c, median_tai), tai_floor)
        for c in codons
    ])

    # raw tau: inverse of tAI (faster decoding = shorter pause)
    raw_tau = 1.0 / tai_values

    # normalize so the mean matches tau_mean
    tau = raw_tau * (tau_mean / np.mean(raw_tau))

    return tau


def shuffled_synonymous_rates(codons, tai_dict, tau_mean=0.06,
                               tai_floor=0.01, seed=None):
    """
    null model: shuffle synonymous codons within each amino acid group,
    then compute tAI-based rates on the shuffled sequence.

    this preserves the amino acid sequence exactly but randomizes
    synonymous codon choice — tests whether specific codon identity
    (not just amino acid identity) matters for translation timing.

    parameters:
      codons    : list of codon strings
      tai_dict  : dict codon -> tAI value
      tau_mean  : desired mean elongation time
      tai_floor : minimum tAI value
      seed      : random seed for reproducibility (None = non-deterministic)

    returns:
      numpy array of elongation times from shuffled codons
    """
    rng = np.random.RandomState(seed) if seed is not None else np.random.RandomState()

    # group positions by amino acid
    # aa_positions: dict aa -> list of position indices
    aa_positions = {}
    for i, codon in enumerate(codons):
        aa = CODON_TABLE.get(codon, None)
        if aa is None or aa == '*':
            # stop codons or unknowns: leave in place
            continue
        if aa not in aa_positions:
            aa_positions[aa] = []
        aa_positions[aa].append(i)

    # build shuffled codon list (start as copy of original)
    shuffled = list(codons)

    # within each AA group, randomly reassign synonymous codons
    for aa, positions in aa_positions.items():
        # collect the codons at these positions
        group_codons = [codons[p] for p in positions]
        # shuffle the codon assignments among these positions
        rng.shuffle(group_codons)
        # write back
        for p, c in zip(positions, group_codons):
            shuffled[p] = c

    # compute tAI rates on the shuffled sequence
    return tai_rates(shuffled, tai_dict, tau_mean=tau_mean, tai_floor=tai_floor)
