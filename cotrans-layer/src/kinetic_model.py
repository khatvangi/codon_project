"""
kinetic model for co-translational folding.

computes contact order, folding rates (plaxco relation),
and time-dependent folding probability during translation.
"""

import numpy as np


def compute_contact_order(contacts, foldon_length):
    """
    compute relative contact order for a foldon segment.

    CO = (1 / (L * N)) * sum(|i - j|)
    where L = foldon_length, N = number of contacts.

    args:
        contacts: list of (i, j) tuples, 0-indexed residue pairs
        foldon_length: number of residues in the foldon segment

    returns:
        float contact order, or None if contacts is empty
    """
    if not contacts:
        return None

    n = len(contacts)
    seq_sep_sum = sum(abs(i - j) for i, j in contacts)
    return seq_sep_sum / (foldon_length * n)


def compute_kf(contact_order):
    """
    estimate folding rate from contact order via plaxco relation.

    ln(kf) = -17.0 * CO + 8.0

    args:
        contact_order: relative contact order (float), or None

    returns:
        kf in s^-1, or None if contact_order is None
    """
    if contact_order is None:
        return None

    return np.exp(-17.0 * contact_order + 8.0)


def cotrans_folding_probability(seg_start, seg_end, kf, ku,
                                 codon_times, tunnel_length=30):
    """
    compute co-translational folding probability at each chain length.

    the foldon segment (seg_start..seg_end) can only begin folding once
    its c-terminal residue clears the ribosome exit tunnel.

    emergence_position = seg_end + tunnel_length

    for each chain length L:
        - if L < emergence_position: P_folded = 0
        - if L >= emergence_position:
            t_avail = sum of codon_times from emergence_position to L (inclusive)
            P_folded = kf/(kf+ku) * (1 - exp(-(kf+ku) * t_avail))

    args:
        seg_start: first residue of foldon (0-indexed)
        seg_end: last residue of foldon (0-indexed)
        kf: folding rate constant (s^-1)
        ku: unfolding rate constant (s^-1)
        codon_times: numpy array, translation time per codon (seconds)
        tunnel_length: ribosome exit tunnel length in residues (default 30)

    returns:
        numpy array of P_folded, same length as codon_times
    """
    n = len(codon_times)
    P = np.zeros(n)

    emergence_position = seg_end + tunnel_length

    # if the foldon never emerges within the transcript, return all zeros
    if emergence_position >= n:
        return P

    # two-state kinetics: rate sum and equilibrium probability
    k_sum = kf + ku
    p_eq = kf / k_sum

    # cumulative time available after emergence.
    # at the emergence position itself, the segment just cleared the tunnel
    # but has had zero time to fold. folding time accumulates from the
    # next codon onward: t_avail(L) = sum(codon_times[emergence_position : L])
    t_avail = 0.0
    for L in range(emergence_position + 1, n):
        t_avail += codon_times[L - 1]
        P[L] = p_eq * (1.0 - np.exp(-k_sum * t_avail))

    return P
