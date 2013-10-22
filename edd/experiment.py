import functools
import collections
import pandas as pa
import numpy as np
from edd import read_bam


class Experiment(object):

    def __init__(self, ip_countd, input_countd, bin_size):
        'should not by instansiated by user-land code'
        self.ipd = ip_countd
        self.inputd = input_countd
        self.bin_size = bin_size

    def aggregate_bins(self, times_bin_size=None, new_bin_size=None):
        if times_bin_size is not None:
            assert times_bin_size > 0
            n = int(times_bin_size)
        elif new_bin_size is not None:
            assert new_bin_size % self.bin_size == 0
            n = int(new_bin_size / self.bin_size)
        else:
            raise Exception("no new bin size given, check api.")
        aipd = read_bam.aggregate_every_n_bins(self.ipd, n)
        ainputd= read_bam.aggregate_every_n_bins(self.inputd, n)
        return Experiment(aipd, ainputd, self.bin_size * n)

    @classmethod
    def load_experiment(cls, chromsizes_path, ip_bam_path, 
            input_bam_path, bin_size=1000, use_multiprocessing=False):
        chromsizes = cls.read_chrom_sizes(chromsizes_path)
        f = functools.partial(read_bam.read_bam_into_bins,
                            chromsizes, bin_size)
        if use_multiprocessing:
            import multiprocessing
            pool = multiprocessing.Pool(processes=2)
            fmap = pool.map
        else:
            fmap = map
        ipd, inputd = fmap(f, [ip_bam_path, input_bam_path])
        return cls(ipd, inputd, bin_size)

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
        ndf = df.copy()
        ndf.input = df.input * input_scale_factor
        return ndf

    def as_data_frame(self, normalize=False):
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

