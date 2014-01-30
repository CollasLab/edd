'''
this module has a somewhat strange flow.
The Experiment class reads bam files and from these raw data
it's really easy to create a dataframe and normalize on those values.

However, other parts of the code (maximum_segments) expects 
an object per bin. Therefore, an utility function coverts 
a df to bin objects. 

'''
import functools
import collections
import itertools
import pandas as pa
import numpy as np
import read_bam, logit
from logbook import Logger
import tempfile
import os
import util
from pybedtools import BedTool
from algorithm.max_segments import GenomeBins, IntervalTest
from algorithm.monte_carlo import MonteCarlo

log = Logger(__name__)

class Experiment(object):
    '''
    classmethod load_experiment reads bam files.
    Instances holds count vectors for IP/INPUT for each chrom
    and knows the bin size (assumed to be fixed)

    The method as_data_frame returns a pandas data frame of
    formatted results. Notice the normalized argument to this
    method.
    '''

    @classmethod
    def load_experiment(cls, chromsizes_path, ip_bam_path, 
            input_bam_path, bin_size=1000, use_multiprocessing=True):
        chromsizes = cls.read_chrom_sizes(chromsizes_path)
        f = functools.partial(read_bam.read_bam_into_bins,
                            chromsizes, bin_size)
        if use_multiprocessing:
            import multiprocessing
            pool = multiprocessing.Pool(processes=2)
            fmap = pool.map
        else:
            fmap = map
        log.notice('loading bam files')
        ipd, inputd = fmap(f, [ip_bam_path, input_bam_path])
        log.notice('done')
        return cls(ipd, inputd, bin_size)


    def __init__(self, ip_countd, input_countd, bin_size):
        'should not by instansiated by user-land code'
        self.ipd = ip_countd
        self.inputd = input_countd
        self.bin_size = bin_size

    def aggregate_bins(self, times_bin_size=None, new_bin_size=None):
        if times_bin_size is not None:
            n = int(times_bin_size)
            assert n > 0
            if n == 1:
                return self
        elif new_bin_size is not None:
            assert new_bin_size % self.bin_size == 0
            n = int(new_bin_size / self.bin_size)
        else:
            raise Exception("no new bin size given, check api.")
        aipd = read_bam.aggregate_every_n_bins(self.ipd, n)
        ainputd= read_bam.aggregate_every_n_bins(self.inputd, n)
        return Experiment(aipd, ainputd, self.bin_size * n)

    def find_smallest_optimal_bin_size(self, nib_lim=0.01, max_ci_diff=0.25):
        for bin_size in itertools.count(1):
            exp = self.aggregate_bins(times_bin_size=bin_size)
            df = exp.as_data_frame()
            ratio_nib = logit.get_nib_ratio(df, max_ci_diff)
            log.notice('testing bin size %d, nib ratio: %.4f' % (bin_size, ratio_nib))
            if ratio_nib <= nib_lim:
                return exp.bin_size
            assert bin_size < 100, "Could not find a suitable bin size."


    @classmethod
    def read_chrom_sizes(cls, chrom_size_filename):
      d = {}
      f = open(chrom_size_filename)
      for line in f:
        chrom, size = line.split()
        if chrom == 'chrom' and size == 'size':
          continue
        d[chrom] = int(size)
      f.close()
      return d

    @classmethod
    def normalize_df(cls, df):
        input_scale_factor = df.ip.sum() / float(df.input.sum())
        log.notice('normalizing input with scale factor: %.2f' % input_scale_factor)
        ndf = df.copy()
        ndf.input = df.input * input_scale_factor
        return ndf

    def as_data_frame(self, normalize=True):
        def chrom_to_df(chrom_name, ip_cnts, input_cnts, bin_size):
            assert len(ip_cnts) == len(input_cnts)
            d = collections.OrderedDict()
            d['chrom'] = chrom_name
            d['start'] = np.arange(len(ip_cnts)) * bin_size
            d['end'] = d['start'] + bin_size
            d['ip'] = ip_cnts
            d['input'] = input_cnts
            return pa.DataFrame(d)
        assert len(self.ipd) == len(self.inputd)
        df = pa.concat([chrom_to_df(c, self.ipd[c], 
            self.inputd[c], self.bin_size)
            for c in self.ipd],
            ignore_index=True)
        if normalize:
            return self.normalize_df(df)
        else:
            return df

    def write_ratios(self, ratio_file):
        log.notice('writing log ratios to %s' % ratio_file)
        df = self.as_data_frame(normalize=True)
        df['ratio'] = np.log(df.ip / df.input).replace(
                [np.inf, -np.inf], np.nan)
        rdf = df.dropna()
        rdf.to_csv(ratio_file, sep='\t', cols=['chrom', 'start', 'end', 'ratio'],
                header=False, index=False)


class BamLoader(object):

    def __init__(self, chrom_size_path, bin_size, neg_score_scale,
                 number_of_processes=4):
        self.chrom_size_path = chrom_size_path
        self.bin_size = bin_size
        self.neg_score_scale = neg_score_scale
        self.bin_size = bin_size
        self.number_of_processes = number_of_processes

    def load_bam(self, ip_name, ctrl_name):
        return Experiment.load_experiment(self.chrom_size_path, ip_name,
                ctrl_name, 1000 if self.bin_size is None else self.bin_size, 
                use_multiprocessing=True)

    def __add_bin_scores(self, r1, r2):
        assert len(r1.index) == len(r2.index)
        assert (r1.index == r2.index).all()
        assert (r1.start == r2.start).all()
        common = r1.copy()
        common.score += r2.score
        return common

    def load_single_experiment(self, ip_name, ctrl_name, gap_file):
        exp = self.load_bam(ip_name, ctrl_name)
        if self.bin_size is None:
            self.bin_size = exp.find_smallest_optimal_bin_size()
            log.notice('Optimal bin size: %d' % self.bin_size)
        else:
            log.notice('Using preset bin size for %s and %s: %d' % (
                ip_name, ctrl_name, self.bin_size))
        odf = exp.aggregate_bins(new_bin_size=self.bin_size).as_data_frame()

        if self.neg_score_scale is None:
            log.notice('Estimating gap penalty')
            self.neg_score_scale = estimate_gap_penalty(odf,
                                                        self.number_of_processes,
                                                        gap_file)
            log.notice('Gap penalty estimated to %.1f' % self.neg_score_scale)
        return logit.ci_for_df(odf, neg_score_scale=self.neg_score_scale)

##############################
# BEGIN ESTIMATE GAP PENALTY #
##############################
import StringIO

def count_stats(xs):
    '''xs is a bedtool instance where the name field holds the bin score'''
    stats = {'DIB': 0, 'EIB': 0}
    for x in xs:
        if float(x.name) > 0:
            stats['EIB'] += 1
        else:
            stats['DIB'] += 1
    return stats

def estimate_gap_penalty(odf, nprocs, gap_file, mc_trials=100, outfile_path=None):
    # gap file marks unalignable regions, something else than
    # gap_penalty TODO : clean this up!
    binscore_df = logit.ci_for_df(odf, neg_score_scale=1)
    bedgraph_path = tempfile.mktemp()
    util.save_bin_score_file(binscore_df, bedgraph_path)
    bg = BedTool(bedgraph_path)
    xs = []
    
    for neg_score_scale in range(2,20):
        log.notice('Testing a gap penalty of %.1f' % neg_score_scale)
        df = logit.ci_for_df(odf, neg_score_scale=neg_score_scale)
        gb = GenomeBins.df_as_bins(df, gap_file)
        max_bin_score = df.score.max()
        observed_result = gb.max_segments(filter_trivial=max_bin_score)
        mc_res = MonteCarlo.run_simulation(gb.chrom_scores, 
                                                         niter=mc_trials, nprocs=nprocs)
        tester = IntervalTest(observed_result, mc_res)
        segments = [segment for (segment, pval) in tester.pvalues() if pval < 0.05]
        peaks_sb = StringIO.StringIO()
        tester.segments_to_bedstream(segments, peaks_sb)
        peaks = BedTool(peaks_sb.getvalue(), from_string=True)
        d = count_stats(bg.intersect(peaks))
        d['gap-penalty'] = neg_score_scale
        xs.append(d)
    genome_wide_stats = count_stats(bg)
    df = pa.DataFrame(xs)
    df['peak_EIB_ratio'] = df.EIB / (df.EIB + df.DIB).astype(float)
    df['global_EIB_coverage'] = df.EIB / float(genome_wide_stats['EIB'])
    df['score'] = df.peak_EIB_ratio**5 * df.global_EIB_coverage
    df.sort('gap-penalty', inplace=True)
    if outfile_path:
        df.to_csv(outfile_path, index=False)
    # TODO perform extra gap_penalty tests for best score +- 0.5
    os.remove(bedgraph_path)
    return df.ix[df.score.argmax()]['gap-penalty']

############################
# END ESTIMATE GAP PENALTY #
############################