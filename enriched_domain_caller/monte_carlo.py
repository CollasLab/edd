import numpy as np
from collections import Counter
from enriched_domain_caller import max_segments
import itertools
import pandas as pa
import sys
from bisect import bisect_left

class MonteCarlo(object):
    """
    Initialized by chromosome sizes (in bins) and number of positive bins.
    A trial performs a random shuffling and a subsequent count of bin sizes.
    """

    def __init__(self, chrom_sizes):
        """
        Arguments:
        - `chrom_sizes`: {<chrom-name>: (num_pos_bins, tot_bins)}
        """
        self._npos, self._ntot = self.__get_gw_stats(chrom_sizes)
        self._chrom_sizes = chrom_sizes
        a = np.empty(self._ntot, dtype=int)
        a[:] = -1
        a[:self._npos] = 1
        self._mc_arr = a

    def __generate_permutation(self):
        a = self._mc_arr#.copy()
        np.random.shuffle(a)
        chrom_lengths = [x[1] for x in self._chrom_sizes.values()]
        # last interval is implicitly known (rest of array), so it's left out
        split_idx = np.cumsum(chrom_lengths[:-1])
        xs = np.split(a, split_idx)
        return dict(zip(self._chrom_sizes, xs))

    def __get_gw_stats(self, chrom_sizes):
        npos = sum(x[0] for x in chrom_sizes.values())
        ntot = sum(x[1] for x in chrom_sizes.values())
        return npos, ntot

    def trial(self):
        segment_value_count = Counter()
        p = self.__generate_permutation()
        for chrom, xs in p.items():
            for segment in max_segments(xs):
                segment_value_count[segment.score] += 1
        return max(segment_value_count.keys())

    def __call__(self, i):
        np.random.seed()
        res = self.trial()
        sys.stdout.write('.')
        sys.stdout.flush()
        return res

def get_sig_limit(obs, mc, fdr_lim):
    '''
    Given two lists of observed/simulated segment score counts,
    Finds the segment score limit that falls below the fdr_limit.
    '''
    def compute_pvalues(obs, mc):
        for o in obs:
            h0_greaterequal = len(mc) - bisect_left(mc, o)
            pvalue = float(h0_greaterequal + 1) / (len(mc) + 1)
            yield pvalue

    def largest_significant_fdr_corr_pval(pvals, fdr_a=.05):
        idx = np.arange(len(pvals)) + 1
        hits = idx * fdr_a / len(pvals) >= pvals
        return np.flatnonzero(hits)[-1]

    obs = np.sort(obs)[::-1]
    mc = np.sort(mc)
    pvals = list(compute_pvalues(obs, mc))
    lim_idx = largest_significant_fdr_corr_pval(pvals, fdr_lim)

    fdr_lim = obs[lim_idx]
    return fdr_lim
