import collections
import itertools
import operator

import logbook
import numpy as np

from enriched_domain_caller import util, gaps, max_segments, score_cutoff

log = logbook.Logger('base')

class GenomeBinScore(object):
    """
    Reads count data.
    """

    def __init__(self, chrom_scores):
        """

        Arguments:
        - `chrom_counts`:
        - `scorefunction`:
        """
        self._chrom_scores = chrom_scores
        self.rev_gaps = None


    @classmethod
    def from_count_file(cls, fname, score_function_name='log2-ratio',
                        normalize=False):
        counts = util.read_counts(fname)
        score_function = util.scoring_functions[score_function_name]
        if normalize:
            input_scale_factor = util.get_input_scale_factor(counts)
            counts = util.normalize_bins(counts, input_scale_factor)
        chromd = collections.defaultdict(list)
        for x in util.score_bins(counts, score_function):
            chromd[x.chrom].append(x)
        return cls(chromd)


    def add_gaps(self, gap_file, drop_gaps_smaller_than):
        """

        Arguments:
        - `gap_file`:
        - `drop_gaps_smaller_than`:
        """
        if gap_file is None:
            return
        assert self.rev_gaps is None
        gs = gaps.read_gap_file(gap_file, drop_gaps_smaller_than)
        d, revdict = gaps.split_on_gaps(self._chrom_scores, gs)
        self._chrom_scores = d
        self.rev_gaps = revdict

    def as_binary(self, min_ratio=0.4):
        self.opt_score = score_cutoff.ScoreCutoff.from_chrom_scores(
            self._chrom_scores).optimize()
        if self.opt_score.ratio > min_ratio:
            log.warn(('Estimated optimal cutoff gives a too high ratio. (%.2f > %.2f'
                     + 'Consider increasing the bin size') % (self.opt_score.ratio, min_ratio))
            self.lim_value = self.opt_score.get_limit_score(min_ratio)
            log.warn('Using non-optimal %.3f as lim value as it gives a ratio of %.2f.' % (self.lim_value, min_ratio))
        else:
            self.lim_value = self.opt_score.get_lim_value_for_ratio(min_ratio)
        score = { True: 1, False: -1 }
        r = {}
        for k, xs in self._chrom_scores.items():
            ys = [util.bed(x.chrom, x.start, x.end,
                           score[x.score > self.lim_value])
                  for x in xs]
            r[k] = ys
        return GenomeBinBinary(r)

class GenomeBinBinary(object):
    """
    """

    def __init__(self, chrom_dict):
        """

        Arguments:
        - `chrom_dict`:
        """
        self._chrom_dict = chrom_dict

    @classmethod
    def bins_to_max_segments_bed(cls, chrom, v):
        xs = max_segments([x.score for x in v])
        return [util.bed(chrom, v[x.from_idx].start,
                    v[x.to_idx].end, x.score) for x in xs]

    def find_maximum_segments(self):
        """
        """
        segments_per_chrom = {}
        for k,v in self._chrom_dict.items():
            segments_per_chrom[k] = self.bins_to_max_segments_bed(k, v)
        return MaximumIntervals(segments_per_chrom)

    def get_stats(self):
        """
        """
        d = {}
        for chrom, spc in self._chrom_dict.items():
            ntot = len(spc)
            npos = len([x.score for x in spc if x.score > 0])
            d[chrom] = (npos, ntot)
        return d

class MaximumIntervals(object):
    """
    """

    def __init__(self, intervals):
        """

        Arguments:
        - `intervals`:
        """
        self._intervals = intervals

    def as_bed(self, output_file, segment_cutoff):
        """

        Arguments:
        - `output_file`:
        - `segment_cutoff`:
        """
        log.notice('Saving significant peaks.')
        cnt = 0
        with open(output_file, 'w') as of:
            for chrom in sorted(self._intervals):
                segments = self._intervals[chrom]
                segments.sort(key=operator.itemgetter(1)) # sort by start index
                for segment in segments:
                    if segment.score >= segment_cutoff:
                        cnt += 1
                        of.write('\t'.join((chrom, str(segment.start),
                                            str(segment.end),
                                            str(segment.score))))
                        of.write('\n')
        log.notice('Done. Wrote %d peaks.' % cnt)

    def filter_trivial(self, lim):
        """
        """
        sumd = sum(len(v) for v in self._intervals.values())
        log.notice('peak filter limit is: %d' % lim)
        log.notice('%8d peaks prior to filtering.' % sumd)

        d = {k:[x for x in xs if x.score > lim]
             for k, xs in self._intervals.items()}

        sumnd = sum(len(v) for v in d.values())
        log.notice('%8d peaks removed after filtering.' % (sumd - sumnd))
        log.notice('%8d peaks remaining after filtering.' % sumnd)
        assert sumd > sumnd
        return MaximumIntervals(d)

    def merge_gaps(self, rev_gaps):
        """

        Arguments:
        - `rev_gaps`:
        """
        if rev_gaps is None:
            # no gaps
            return self
        merged_intervals = gaps.join_gaps(self._intervals, rev_gaps)
        return MaximumIntervals(merged_intervals)

    def get_interval_scores(self):
        all_intervals = itertools.chain.from_iterable(self._intervals.values())
        return np.array([x.score for x in all_intervals])
