from algorithm.max_segments import GenomeBins, IntervalTest
from algorithm.monte_carlo import MonteCarlo
from pybedtools import BedTool
import StringIO
import logit
import tempfile
import util
import os
from math import sqrt
from logbook import Logger
import itertools
import pandas as pa
import scipy.stats
log = Logger(__name__)

def corrcoeff(odf):
    df = pa.DataFrame(dict(left=odf.score[:-1].values,
                           right=odf.score[1:].values)).dropna()
    return scipy.stats.spearmanr(df.left, df.right)[0]
    
def bin_size(orig_exp, nib_lim=0.01, max_ci_diff=0.25, min_corcoef=0.3):
    for bin_size in itertools.count(1):
        exp = orig_exp.aggregate_bins(times_bin_size=bin_size)
        df = logit.ci_for_df(exp.as_data_frame(), max_ci_diff,
                             neg_score_scale=1, extrapolate_low_info_bins=False)
        ratio_nib = logit.get_nib_ratio(df)
        pvar = corrcoeff(df)
        log.notice('testing bin size %d, nib ratio: %.4f, spearmanr: %.3f' % (bin_size, ratio_nib, pvar))
        if ratio_nib <= nib_lim and pvar > min_corcoef:
            return exp.bin_size
        assert bin_size < 100, "Could not find a suitable bin size."

def golden_section_search(f, left, mid, right, precision):
    phi = (1 + sqrt(5))/2
    resphi = 2 - phi

    def g(l, m, r):
        if abs(l - r) < precision:
            return (l + r)/2.0
        # Create a new possible center, in the area between c and b, pushed against c
        mr = m + resphi*(r - m)
        log.notice('searching [l=%.2f, m=%.2f, mr=%.2f, r=%.2f]' % (l, m, mr, r))
        mr_res = f(mr)
        m_res = f(m)

        if f(mr) > f(m):
            log.notice('f(mr)=%.2f > f(m)=%.2f -> g(m, mr, r)' % (mr_res, m_res))
            return g(m, mr, r)
        else:
            log.notice('f(mr)=%.2f <= f(m)=%.2f -> g(mr, m, l)' % (mr_res, m_res))
            return g(mr, m, l)
    return g(left, mid, right)

        
class GapPenalty(object):



    def __init__(self, orig_bins, bedgraph_path, nprocs, gap_file,
                 mc_trials, pval_lim):
        self.orig_bins = orig_bins
        self.bedgraph_path = bedgraph_path
        self.bins_bedtool = BedTool(self.bedgraph_path)
        self.nprocs = nprocs
        self.gap_file = gap_file
        self.mc_trials = mc_trials
        self.pval_lim = pval_lim
        self.__cache = {}
        self.genome_wide_stats = self.count_stats(self.bins_bedtool)

    def cleanup(self):
        os.remove(self.bedgraph_path) # ugly I know

    @classmethod
    def count_stats(self, xs):
        '''xs is a bedtool instance where the name field holds the bin score'''
        stats = {'DIB': 0, 'EIB': 0}
        for x in xs:
            if float(x.name) > 0:
                stats['EIB'] += 1
            else:
                stats['DIB'] += 1
        return stats
        
    @classmethod
    def instantiate(cls, odf, nprocs, gap_file, mc_trials, pval_lim):
        binscore_df = logit.ci_for_df(odf, neg_score_scale=1)
        binscore_gb = GenomeBins.df_as_bins(binscore_df, gap_file)
        bedgraph_path = tempfile.mktemp()
        util.save_bin_score_file(binscore_df, bedgraph_path)
        rval = cls(binscore_gb, bedgraph_path, nprocs, gap_file, mc_trials, pval_lim)
        return rval


    def search(self, left=2.0, mid=10.0, right=30, precision=0.1):
        return golden_section_search(self.comp_score, left, mid, right,
                                     precision)
        
    def comp_score(self, gap_penalty):
        '''compute_score_given_gap_penalty'''
        if gap_penalty in self.__cache:
            return self.__cache[gap_penalty]['score']
        
        gb = self.orig_bins.scale_neg_scores(gap_penalty)
        observed_result = gb.max_segments()
        mc_res = MonteCarlo.run_simulation(gb.chrom_scores, 
                                           niter=self.mc_trials, nprocs=self.nprocs)
        tester = IntervalTest(observed_result, mc_res)
        segments = [segment for (segment, pval) in tester.pvalues()
                    if pval < self.pval_lim]
        # TODO use bx.python instead of pybedtools
        peaks_sb = StringIO.StringIO()
        tester.segments_to_bedstream(segments, peaks_sb)
        peaks = BedTool(peaks_sb.getvalue(), from_string=True)
        d = self.count_stats(self.bins_bedtool.intersect(peaks))
        d['gap-penalty'] = gap_penalty
        d['peak_EIB_ratio'] = d['EIB'] / float(d['EIB'] + d['DIB'])
        d['global_EIB_coverage'] = d['EIB'] / float(self.genome_wide_stats['EIB'])
        d['score'] = d['peak_EIB_ratio']**5 * d['global_EIB_coverage']
        log.notice('Gap penalty of %.2f gives a score of %.3f' % (gap_penalty, d['score']))
        self.__cache[gap_penalty] = d
        return d['score']
