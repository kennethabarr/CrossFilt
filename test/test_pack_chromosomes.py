"""
Unit tests for liftover_functions.pack_chromosomes.

Run with:
    pytest test/test_pack_chromosomes.py
"""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from crossfilt.lib.liftover_functions import pack_chromosomes


class _Stat:
    """Minimal stand-in for pysam index-statistics objects."""
    def __init__(self, contig, total):
        self.contig = contig
        self.total  = total


# ------------------------------------------------------------------ #
#  Empty / degenerate inputs                                          #
# ------------------------------------------------------------------ #

def test_empty_stats_returns_empty():
    assert pack_chromosomes([], n_files=4) == []


def test_all_zero_counts_returns_empty():
    stats = [_Stat('chr1', 0), _Stat('chr2', 0)]
    assert pack_chromosomes(stats, n_files=4) == []


# ------------------------------------------------------------------ #
#  Bin count capping                                                  #
# ------------------------------------------------------------------ #

def test_fewer_chroms_than_bins():
    """Never return more bins than there are non-empty chromosomes."""
    stats = [_Stat('chr1', 100), _Stat('chr2', 50)]
    result = pack_chromosomes(stats, n_files=10)
    assert len(result) == 2
    all_chroms = [c for bin_ in result for c in bin_]
    assert sorted(all_chroms) == ['chr1', 'chr2']


def test_single_chrom_single_bin():
    stats = [_Stat('chr1', 1000)]
    result = pack_chromosomes(stats, n_files=1)
    assert result == [['chr1']]


def test_exact_one_chrom_per_bin():
    stats = [_Stat(f'chr{i}', (i + 1) * 100) for i in range(4)]
    result = pack_chromosomes(stats, n_files=4)
    assert len(result) == 4
    all_chroms = sorted(c for bin_ in result for c in bin_)
    assert all_chroms == sorted(s.contig for s in stats)
    for bin_ in result:
        assert len(bin_) == 1


# ------------------------------------------------------------------ #
#  All chromosomes are distributed (no reads dropped)                #
# ------------------------------------------------------------------ #

def test_all_chroms_assigned():
    stats = [_Stat(f'chr{i}', (i + 1) * 500) for i in range(10)]
    result = pack_chromosomes(stats, n_files=4)
    all_chroms = sorted(c for bin_ in result for c in bin_)
    expected   = sorted(s.contig for s in stats)
    assert all_chroms == expected


def test_no_duplicate_assignments():
    stats = [_Stat(f'chr{i}', i * 100 + 1) for i in range(8)]
    result = pack_chromosomes(stats, n_files=3)
    all_chroms = [c for bin_ in result for c in bin_]
    assert len(all_chroms) == len(set(all_chroms))


# ------------------------------------------------------------------ #
#  Balance (LPT greedy should produce near-optimal balance)          #
# ------------------------------------------------------------------ #

def test_balanced_equal_chroms():
    """Four equal chromosomes into two bins should be perfectly balanced."""
    stats = [_Stat(f'chr{i}', 1000) for i in range(4)]
    result = pack_chromosomes(stats, n_files=2)
    assert len(result) == 2
    totals = [sum(1000 for _ in bin_) for bin_ in result]
    assert totals[0] == totals[1] == 2000


def test_lpt_handles_skewed_distribution():
    """
    One huge chromosome and many small ones.
    LPT assigns the largest first, so chrBig lands alone in its own bin.
    The remaining bins share the small chromosomes and should be balanced.
    """
    stats = [_Stat('chrBig', 10000)] + [_Stat(f'chr{i}', 100) for i in range(20)]
    result = pack_chromosomes(stats, n_files=4)
    assert len(result) == 4
    # chrBig must be alone (LPT assigns it first to an empty bin)
    big_bins = [bin_ for bin_ in result if 'chrBig' in bin_]
    assert len(big_bins) == 1
    assert big_bins[0] == ['chrBig']
    # the remaining three bins share 20 * 100 = 2000 reads; each should
    # receive between 5 and 8 chromosomes (greedy near-equal split)
    other_bins = [bin_ for bin_ in result if 'chrBig' not in bin_]
    for bin_ in other_bins:
        assert 5 <= len(bin_) <= 8


def test_known_assignment():
    """
    Hand-verified LPT assignment for a small case.
    chroms: A=5000, B=4000, C=3000, D=2000  into 2 bins
    LPT: assign A(5000) → bin0, B(4000) → bin1, C(3000) → bin0 (5000<4000? no)
         actually: bin0=5000, assign B→bin1 (4000), assign C→bin0 (5000→8000 no)
         C→bin1 (4000→7000), D→bin0 (5000→7000)
    So both bins end up at 7000 each.
    """
    stats = [_Stat('A', 5000), _Stat('B', 4000), _Stat('C', 3000), _Stat('D', 2000)]
    result = pack_chromosomes(stats, n_files=2)
    totals = sorted(
        sum(next(s.total for s in stats if s.contig == c) for c in bin_)
        for bin_ in result
    )
    assert totals == [7000, 7000]
