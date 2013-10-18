from libc.stdlib cimport calloc, free
from libc.string cimport strncpy
import numpy as np
cimport numpy as np
cimport pysam.csamtools as csam

def read_bam_into_bins(chrom_sizes, bin_size, bam_filename):
  b = BamCounter(chrom_sizes, bam_filename, bin_size)
  b.process_bam()
  return b.chrom_bins

cdef class BamCounter:
  cdef:
    object chrom_sizes
    public object chrom_bins
    double **cb
    size_t cb_len 
    csam.Samfile fp 
    size_t bin_size

  def __cinit__(self, chrom_sizes, bam_filename, bin_size):
    self.fp = csam.Samfile(bam_filename)
    self.cb_len = self.fp.nreferences
    self.cb = <double**> calloc(self.cb_len, sizeof(double*));

  def __init__(self, chrom_sizes, bam_filename, bin_size):
    self.chrom_sizes = chrom_sizes
    self.bin_size = bin_size
    self.chrom_bins = {}
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] a

    # add array bin pointers to fast lookup array
    for i in range(self.cb_len):
      chrom_name = self.fp.getrname(i) 
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
      csam.bam1_t *b
      size_t start, cov, sidx
      size_t nreads = 0
      size_t start_to_bin_end
      double r 
    while self.fp.cnext() >= 0:
      b = self.fp.getCurrent()
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
    return nreads;

  def __dealloc__(self):
    free(self.cb);
