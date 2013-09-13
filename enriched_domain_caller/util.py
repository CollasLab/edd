import matplotlib.pylab as plt
import os
import numpy as np
from collections import namedtuple
import math
import operator
import sys
from logbook import Logger
#import ipdb
verbose_dir = None
log = Logger('base')
hg19_chromfilter = set(['chrY', 'chrX', 'chr13', 'chr12', 'chr11', 'chr10', 'chr17', 'chr16', 'chr15', 'chr14', 'chr19', 'chr18', 'chr22', 'chr20', 'chr21', 'chr7', 'chr6', 'chr5', 'chr4', 'chr3', 'chr2', 'chr1', 'chr9', 'chr8'])
bed = namedtuple('BedGraph', 'chrom start end score')
bin = namedtuple('Bin', 'chrom start end ip_count input_count')

def log2_score(pos, neg):
    if pos == 0 or neg == 0:
        return -5 # default negative for missing values
    return math.log(float(pos) / neg)

def ci_lower_bound(pos, neg):
    '''
    computes lower bound of 95% confidence interval for true binomial proportion.

    sources:
    http://www.evanmiller.org/how-not-to-sort-by-average-rating.html
    http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    '''
    n = pos + neg
    if n == 0:
        return 0
    z = 1.96 # for a 95% confidence interval
    pratio = float(pos) / n
    return (pratio + z*z/(2*n) - z * math.sqrt((pratio*(1-pratio)+z*z/(4*n))/n))/(1+z*z/n)

scoring_functions = {
    'log2-ratio': log2_score,
    'CI-lower': ci_lower_bound
    }

def get_input_scale_factor(bins):
    '''
    TODO: move somewhere
    '''
    ip_sum = input_sum = 0
    for x in bins:
        ip_sum += x.ip_count
        input_sum += x.input_count
    scale_fac = float(ip_sum) / input_sum
    log.notice('input scale factor is: %.2f' % scale_fac)
    return scale_fac

def normalize_bins(bins, input_scale_factor):
    '''
    TODO: move somewhere
    '''
    xs = []
    for y in bins:
        x = bin(y.chrom, y.start, y.end, y.ip_count,
                y.input_count * input_scale_factor)
        xs.append(x)
    return xs

def score_bins(bins, scorefunc):
    bs = []
    for x in bins:
        score = scorefunc(x.ip_count, x.input_count)
        bs.append(bed(x.chrom, x.start, x.end, score))
    return bs

def read_counts(fname):
    bins = []
    with open(fname) as f:
        for line in f:
            parts = line.split()
            start, end = int(parts[1]), int(parts[2])
            chrom = parts[0]
            ip_cnt = float(parts[3])
            input_cnt = float(parts[4])
            if ip_cnt == 0 or input_cnt == 0:
                # allmost all bins have some level of background.
                # remove those without
                ip_cnt = input_cnt = 0
            bins.append(bin(chrom, start, end, ip_cnt, input_cnt))
    assert len(bins) > 0
    return bins
