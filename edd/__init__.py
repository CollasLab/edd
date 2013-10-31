import pysam as __pysam
import max_segments
import read_bam

import experiment

load_experiment = experiment.Experiment.load_experiment 
class NoPeaksException(RuntimeError):
    pass

def parse_bin_size_as_single_number(s):
    if len(s) > 2 and s[-2:].lower() == 'kb':
        return int(s[:-2]) * 1000
    else:
        return int(s)

