import StringIO
import pysam as __pysam
import collections
import read_bam
import experiment
from estimate_gap_penalty import estimate_gap_penalty
import util
from logbook import Logger
import algorithm
log = Logger(__name__)

load_experiment = experiment.Experiment.load_experiment


def df_as_bins(df, gap_file):
    '''
    converts a already scored df to an object
    containing a dict of bins per chromosome, separated by gaps.
    '''
    chromd = collections.defaultdict(list)
    for _, x in df.iterrows():
        b = util.bed(x['chrom'], x['start'], x['end'],
                x['score'])
        chromd[b.chrom].append(b)
    return algorithm.GenomeBins.with_gaps(chromd, gap_file)

class BamLoader(object):

    def __init__(self, chrom_size_path, bin_size, neg_score_scale,
                 number_of_processes=4):
        self.chrom_size_path = chrom_size_path
        self.bin_size = bin_size
        self.neg_score_scale = neg_score_scale
        self.bin_size = bin_size
        self.number_of_processes = number_of_processes

    def load_bam(self, ip_name, ctrl_name):
        return edd.load_experiment(self.chrom_size_path, ip_name,
                ctrl_name, 1000 if self.bin_size is None else self.bin_size, 
                use_multiprocessing=True)

    def __add_bin_scores(self, r1, r2):
        assert len(r1.index) == len(r2.index)
        assert (r1.index == r2.index).all()
        assert (r1.start == r2.start).all()
        common = r1.copy()
        common.score += r2.score
        return common

    def load_single_experiment(self, ip_name, ctrl_name):
        exp = self.load_bam(ip_name, ctrl_name)
        if self.bin_size is None:
            self.bin_size = exp.find_smallest_optimal_bin_size()
            log.notice('Optimal bin size: %d' % self.bin_size)
        else:
            log.notice('Using preset bin size for %s and %s: %d' % (
                ip_name, ctrl_name, self.bin_size))
        odf = exp.aggregate_bins(new_bin_size=self.bin_size).as_data_frame()
        if self.neg_score_scale is None:
            selg.neg_score_scale = estimate_gap_penalty(odf,
                                                        number_of_processes)
        return edd.logit.ci_for_df(odf, neg_score_scale=self.neg_score_scale)

