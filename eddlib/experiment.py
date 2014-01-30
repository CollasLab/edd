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
from edd import read_bam, logit
from logbook import Logger
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


