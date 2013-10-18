import functools
import multiprocessing
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
    f = functools.partials(read_bam.read_bam_into_bins,
            chromsizes, bin_size)
    ipd, inputd = pool.map(f, [ip_bam_path, input_bam_path])
    return dict(ip=ipd, input=inputd)
