from libc.stdlib cimport calloc, free
from libc.string cimport strncpy
import numpy as np
cimport numpy as np
from pysam.csamtools cimport samfile_t, bam1_t, samopen, bam_init1, bam_destroy1, samread

DEF MAX_NAME_LEN = 50

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

cdef class BamCounter:
  cdef:
    object chrom_sizes
    public object chrom_bins
    double **cb
    size_t cb_len 
    samfile_t *fp
    size_t bin_size

  def __cinit__(self, chrom_sizes, bam_filename, bin_size):
    self.fp = samopen(bam_filename, "rb", NULL)
    if self.fp == NULL:
      raise IOError("Failed to open BAM file %s\n" % bam_filename);
    self.cb_len = self.fp.header.n_targets
    self.cb = <double**> calloc(self.cb_len, sizeof(double*));

  def __init__(self, chrom_sizes, bam_filename, bin_size):
    self.chrom_sizes = chrom_sizes
    self.bin_size = bin_size
    self.chrom_bins = {}
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] a

    # add array bin pointers to fast lookup array
    for i in range(self.fp.header.n_targets):
      chrom_name = self.fp.header.target_name[i] 
      if not chrom_name in self.chrom_sizes:
        self.cb[i] = NULL
      else:
        chrom_size = self.chrom_sizes[chrom_name]
        nbins = int(chrom_size / bin_size)
        if chrom_size % bin_size != 0:
          nbins += 1
        a = np.zeros(nbins, dtype=np.float64)
        self.chrom_bins[chrom_name] = a
        self.cb[i] = <double *> a.data

  def process_bam(self):
    cdef:
      bam1_t *b = bam_init1()
      size_t start, cov, sidx
      size_t nreads = 0
      size_t start_to_bin_end
      double r 
    while (samread(self.fp, b) >= 0):
      if (b.core.flag & 0x4 or 
          self.cb[b.core.tid] == NULL or 
          b.core.pos == 0):
        continue
      nreads += 1
      start = b.core.pos - 1 # bam 1-based, bed is 0-based
      sidx = start / self.bin_size
      start_to_bin_end = self.bin_size - (start % self.bin_size)
      if b.core.l_qseq < start_to_bin_end:
        # the whole sequence fits within a bin
        self.cb[b.core.tid][sidx] += 1
      else:
        # the sequence spans two bins
        r = start_to_bin_end / <double>b.core.l_qseq;
        self.cb[b.core.tid][sidx] += r
        self.cb[b.core.tid][sidx + 1] += (1 - r)
    bam_destroy1(b);
    return nreads;

  def __dealloc__(self):
    free(self.cb);
