import numpy as np
from collections import Counter
from enriched_domain_caller import max_segments

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
        segment_value_count = Counter()
        p = self.__generate_permutation()
        for chrom, xs in p.items():
            for segment in max_segments(xs):
                segment_value_count[segment.score] += 1
        return segment_value_count

    def __call__(self, i):
        print 'iteration:', i
        return self.trial()
