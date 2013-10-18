import pysam as __pysam
import max_segments
import read_bam

class NoPeaksException(RuntimeError):
    pass
