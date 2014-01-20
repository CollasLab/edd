import operator
import numpy as np
import collections
import toolz
from util import Bedgraph
from chrom_max_segments import max_segments
from edd.algorithm import gaps
import logbook
from rpy2.robjects.packages import importr

stats = importr('stats')
log = logbook.Logger(__name__)

class GenomeBins(object):
    '''
    Genome wide bins with scores.
    '''

    def __init__(self, chrom_bins):
        '''
        Keyword arguments:
        chrom_bins -- dict of names->bins where a name may span a 
                        complete chromosome or just a part
        '''
        self.chrom_bins = chrom_bins
        self.chrom_scores = {k:np.array([x.score for x in v], dtype=float) 
                for k, v in chrom_bins.items()}

    @classmethod
    def with_gaps(cls, chrom_bins, gap_file):
        if gap_file is not None:
            g = gaps.read_gap_file(gap_file)
            chrom_bins, _ = gaps.split_on_gaps(chrom_bins, g)
        return cls(chrom_bins)

    def max_segments(self, filter_trivial=0):
        def chrom_max(chrom):
            xs = max_segments(self.chrom_scores[chrom])
            bins = self.chrom_bins[chrom]
            return [Bedgraph(bins[x.from_idx].chrom, bins[x.from_idx].start,
                        bins[x.to_idx].end, x.score) for x in xs]
        segments_per_chrom = collections.defaultdict(list)
        for x in toolz.concat(chrom_max(k) for k in self.chrom_bins):
            if x.score > filter_trivial:
                segments_per_chrom[x.chrom].append(x)
        num_intervals = 0
        for xs in segments_per_chrom.values():
            xs.sort(key=operator.attrgetter('start'))
            num_intervals += len(xs)
        log.notice('Removed trivial intervals with score less than %.4f.' % filter_trivial)
        log.notice('%d intervals (potential peaks) remaining.' % num_intervals)
        return MaxIntervals(segments_per_chrom)

class MaxIntervals(object):

    def __init__(self, intervals):
        self._intervals = intervals

    def as_bed(self, output_file, segment_cutoff):
        """
        Arguments:
        - `output_file`:
        - `segment_cutoff`: saving peaks with a score higher than cutoff
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

    def all_scores(self, filter_trivial=0):
        return np.array([x.score
            for x in toolz.concat(self._intervals.itervalues())
            if x.score > filter_trivial])

class IntervalTest(object):

    def __init__(self, max_intervals, mc_res):
        self.max_intervals = max_intervals
        self.mc_res = np.sort(mc_res)
        self._pvalues = None
        self._pvalues = None

    def pvalues(self):
        res = []
        for xs in self.max_intervals._intervals.values():
            pos = len(self.mc_res) - np.searchsorted(self.mc_res,
                    [x.score for x in xs])
            pvals = (pos + 1) / float(len(self.mc_res) + 1)
            for x, pval in zip(xs, pvals):
                res.append((x, pval))
        self._pvalues = res
        return res

    def qvalues(self, below=0.1):
        if not self._pvalues:
            self.pvalues()
        pvals = [x[1] for x in self._pvalues]
        qvals = list(stats.p_adjust(pvals, method='fdr'))
        res = [(q,p,x) for (q, (x,p)) in zip(qvals, self._pvalues)
            if q < below]
        self._qvalues = res
        log.notice('got %d peaks with qvalue below %.2f. From %d possible.' % (
            len(res), below, len(pvals)))
        return res

    def as_bed(self, output_file):
        if self._qvalues is None:
            self.qvalues()
        log.notice('Saving significant peaks.')
        with open(output_file, 'w') as of:
            segments = [x for (_, _, x) in self._qvalues]
            segments.sort(key=operator.itemgetter(0,1)) # sort by chrom, then start
            for segment in segments:
                of.write('\t'.join((segment.chrom, str(segment.start),
                                    str(segment.end),
                                    str(segment.score))))
                of.write('\n')
        log.notice('Done. Wrote %d peaks.' % len(self._qvalues))

