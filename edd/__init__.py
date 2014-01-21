import pysam as __pysam
import collections
import read_bam
import experiment
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

