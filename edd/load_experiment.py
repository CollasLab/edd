import functools
import multiprocessing
import collections
import pandas as pa
import numpy as np
from edd import read_bam

def read_chrom_sizes(chrom_size_filename):
  d = {}
  f = open(chrom_size_filename)
  for line in f:
    chrom, size = line.split()
    if chrom == 'chrom' and size == 'size':
      continue
    d[chrom] = int(size)
  f.close()
  return d

def load_experiment(chromsizes_path, ip_bam_path, input_bam_path, bin_size=1000):
    pool = multiprocessing.Pool(processes=2)
    chromsizes = read_chrom_sizes(chromsizes_path)
    f = functools.partial(read_bam.read_bam_into_bins,
            chromsizes, bin_size)
    ipd, inputd = pool.map(f, [ip_bam_path, input_bam_path])
    return ipd, inputd

def experiment_as_df(ipd, inputd, bin_size):
    ''' `d` is a dict of {'ip', 'input'}
    each variable contains a dict of chrom and bin scores.

    purpose: use this + bin_size to get a df with
    chrom start end ip_cnt input_cnt
    '''
    assert len(ipd) == len(inputd)
    return pa.concat([chrom_to_df(c, ipd[c], inputd[c], bin_size)
        for c in ipd])

def chrom_to_df(chrom_name, ip_cnts, input_cnts, bin_size):
    assert len(ip_cnts) == len(input_cnts)
    d = collections.OrderedDict()
    d['chrom'] = chrom_name
    d['start'] = np.arange(len(ip_cnts)) * bin_size
    d['end'] = d['start'] + bin_size
    d['ip'] = ip_cnts
    d['input'] = input_cnts
    return pa.DataFrame(d)
