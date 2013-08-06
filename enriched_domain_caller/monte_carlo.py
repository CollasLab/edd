import numpy as np
from collections import defaultdict
from enriched_domain_caller import max_segments
import itertools
import pandas as pa
import sys

def create_manual_counter(xs=None):
    c = defaultdict(int)
    if xs is not None:
        for x in xs:
            c[x] += 1
    return c

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
        np.random.shuffle(self._mc_arr)
        chrom_lengths = [x[1] for x in self._chrom_sizes.values()]
        # last interval is implicitly known (rest of array), so it's left out
        split_idx = np.cumsum(chrom_lengths[:-1])
        xs = np.split(self._mc_arr, split_idx)
        return dict(zip(self._chrom_sizes, xs))

    def __get_gw_stats(self, chrom_sizes):
        npos = sum(x[0] for x in chrom_sizes.values())
        ntot = sum(x[1] for x in chrom_sizes.values())
        return npos, ntot

    def trial(self):
        segment_value_count = create_manual_counter()
        p = self.__generate_permutation()
        for chrom, xs in p.items():
            for segment in max_segments(xs):
                segment_value_count[segment.score] += 1
        return segment_value_count

    def __call__(self, i):
        res = self.trial()
        sys.stdout.write('.')
        sys.stdout.flush()
        return res

def get_sig_limit(obs, mc, fdr_lim):
    '''
    Given two count dicts of observed/simulated segment score counts,
    Finds the segment score limit that falls below the fdr_limit.
    FDR(score) is computed as sum(freq(H0_i)) / sum(freq(Obs_i))
      for i from score to max_score
    '''
    print 'Filtering significant peaks.',
    def cumulative_prct(xs):
        return xs[::-1].cumsum()[::-1] / float(xs.sum())

    def make_dense(xs):
        i = xs.index.values
        reps = np.concatenate((i[1:] - i[:-1], [1]))
        xs = xs.repeat(reps).values
        idx = np.arange(i.min(), i.max() + 1)
        return pa.Series(xs, index=idx)

    def increase_length(a, l):
        assert len(a) < l
        start = len(a) + 1
        end = l
        idx = np.arange(start, end + 1)
        s = pa.Series(np.zeros(end - start + 1), index=idx)
        return pa.concat([a, s])
    scores = [x.score for x in itertools.chain.from_iterable(obs.values())]
    h0 = pa.Series(mc)
    s = pa.Series(create_manual_counter(scores))
    obs_prct = make_dense(cumulative_prct(s))
    h0_prct = make_dense(cumulative_prct(h0))
    if len(h0_prct) > len(obs_prct):
        obs_prct = increase_length(obs_prct, len(h0_prct))
    else:
        h0_prct = increase_length(h0_prct, len(obs_prct))

    fdr = h0_prct / obs_prct
    lim = fdr[(fdr < fdr_lim)].index[0]
    print 'Done'
    return lim
