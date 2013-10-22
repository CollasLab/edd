import pysam as __pysam
import max_segments
import read_bam

import experiment

load_experiment = experiment.Experiment.load_experiment 
class NoPeaksException(RuntimeError):
    pass
