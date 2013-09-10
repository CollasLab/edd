import itertools
import matplotlib.pylab as plt
import os
import numpy as np
from max_segments import max_segments
from collections import namedtuple, defaultdict
import math
import operator
import multiprocessing
import monte_carlo as mcarlo
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

def monte_carlo(chrom_sizes, niter=4, nprocs=4):
    mc = mcarlo.MonteCarlo(chrom_sizes)
    sys.stdout.write('Performing %d monte carlo trials: ' % niter)
    sys.stdout.flush()
    if nprocs > 1:
        m = multiprocessing.Pool(nprocs)
        xs = m.map(mc, range(niter))

    else:
        xs = [mc(i) for i in range(niter)]
    sys.stdout.write('\nDone\n')
    return np.sort(xs)

def score_file_to_bins(fname, negative_score_factor=1.0):
    '''ipython'''
    d = defaultdict(list)
    for line in open(fname):
        parts = line.split()
        score = float(parts[3])
        if score < 0:
            score *= negative_score_factor
        x = bed(parts[0], int(parts[1]), int(parts[2]), score)
        d[x.chrom].append(x)
    return d


def bins_to_max_segments_bed(chrom, v):
    xs = max_segments([x.score for x in v])
    return [bed(chrom, v[x.from_idx].start,
                v[x.to_idx].end, x.score) for x in xs]

def obs_results(scores_per_chrom):
    segments_per_chrom = {}
    for k,v in scores_per_chrom.items():
        segments_per_chrom[k] = bins_to_max_segments_bed(k, v)
    return segments_per_chrom

def filter_smaller_than_Nsd_from_mean(d, N):
    scores = np.array([x.score for x in itertools.chain.from_iterable(d.values())])
    log.notice('%d potential peaks are FDR corrected (%.2f STD[%.2f] from mean[%.2f])' % (
        len(scores), N, scores.std(), scores.mean()))
    assert scores.std() > 1.
    lim = scores.mean() + N * scores.std()
    return filter_smaller_than_lim(d, lim)

def filter_smaller_than_lim(d, lim):
    sumd = sum(len(v) for v in d.values())
    log.notice('peak filter limit is: %d' % lim)
    log.notice('%8d peaks prior to filtering.' % sumd)
    nd = {k:[x for x in v if x.score > lim]
          for k, v in d.items()}
    sumnd = sum(len(v) for v in nd.values())
    log.notice('%8d peaks removed after filtering.' % (sumd - sumnd))
    log.notice('%8d peaks remaining after filtering.' % sumnd)
    assert sumd > sumnd
    return nd

def read_counts(bedgraph, legal_chroms):
    bins = []
    for line in bedgraph:
        parts = line.split()
        start, end = int(parts[1]), int(parts[2])
        chrom = parts[0]
        if legal_chroms is None or chrom in legal_chroms:
            ip_cnt = float(parts[3])
            input_cnt = float(parts[4])
            if ip_cnt == 0 or input_cnt == 0:
                # allmost all bins have some level of background.
                # remove those without
                ip_cnt = input_cnt = 0
            bins.append(bin(chrom, start, end, ip_cnt, input_cnt))
    assert len(bins) > 0
    return bins

def get_input_scale_factor(bins):
    ip_sum = input_sum = 0
    for x in bins:
        ip_sum += x.ip_count
        input_sum += x.input_count
    scale_fac = float(ip_sum) / input_sum
    log.notice('input scale factor is: %.2f' % scale_fac)
    return scale_fac

def normalize_bins(bins, input_scale_factor):
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

def compute_neg_score(sbins, max_pos_ratio):
    p, n = 0, 0
    for x in sbins:
        if x.score > 0:
            p += 1
        else:
            n += 1
    cur_ratio = float(p) / (p + n)
    if cur_ratio <= max_pos_ratio:
        log.info('cur_ratio <= max_pos_ratio --> -1')
        neg_score = -1
    else:
        r = 1 - max_pos_ratio
        fac = (p*r) / (n*(1-r))
        neg_score = -1 * fac
    assert neg_score <= -1
    log.notice('neg score factor is: %.2f' % neg_score)
    return neg_score

def get_limit_score(bs, pos_bin_score_ratio):
    bs = list(bs)
    bs.sort(reverse=True, key=operator.attrgetter('score'))
    lim_score_idx = int(len(bs) * pos_bin_score_ratio)
    while True:
        lim_score = bs[lim_score_idx].score
        if lim_score > 0:
            break
        lim_score_idx -= 1
    return lim_score

def get_bedgraph_list(bs, lim_score):
    bs = []
    for x in bins:
        val = 1 if x.score >= lim_score else -1
        b = bed(x.chrom, x.start, x.end, val)
        bs.append(b)
    return b

def information_score(bins):
    '''
    returns a score that tries to say something
    about coverage to quality.
    (we want to cover many high quality bins)
    '''
    npos = bins.sum()
    r = float(npos) / len(bins)
    expected = r**2 * (len(bins) - 1)
    observed = np.logical_and(bins[:-1], bins[1:]).sum()
    information_content = math.log(observed / float(expected))
    return information_content * npos


def optimize_score_cutoff(scores):
    '''
    we require the cutoff to be between 0 and mean pos score
    '''
    lim = scores[scores > 0].mean()
    log.notice('searching for optimal pos bin ratio lim between 0 and %.3f.' % lim)
    xs = np.linspace(0, lim, 70)
    ys = np.array([information_score(scores > pos_cutoff)
                   for pos_cutoff in xs])
    if verbose_dir is not None:
        opath = os.path.join(verbose_dir, 'information_content.png')
        plt.clf()
        plt.plot(xs, ys)
        plt.title('log2(obs / expected) * npos')
        plt.savefig(opath)
    lim_value = xs[ys.argmax()]
    orig_pos_ratio = (scores > 0).sum() / float(len(scores))
    pos_ratio = (scores > lim_value).sum() / float(len(scores))
    log.notice('Original positive bin ratio is %.2f' % orig_pos_ratio)
    log.notice('Adjusted positive bin ratio is %.2f' % pos_ratio)
    return lim_value

def read_scores(bedgraph, legal_chroms, scorefunc):
    if isinstance(scorefunc, str):
        scorefunc = scoring_functions[scorefunc]
    bins = read_counts(bedgraph, legal_chroms)
    assert len(bins) > 0
    input_scale_factor = get_input_scale_factor(bins)
    nbins = normalize_bins(bins, input_scale_factor)
    sbins = score_bins(nbins, scorefunc)
    lim_score = optimize_score_cutoff(np.array([x.score for x in sbins]))
    pos_score = 1
    neg_score = -1 #compute_neg_score(sbins, max_pos_ratio=pos_bin_ratio)

    chromd = defaultdict(list)
    for x in sbins:
        bin_score = 1 if x.score > lim_score else neg_score
        chromd[x.chrom].append(bed(x.chrom, x.start, x.end, bin_score))
    return chromd

def write_segments(of, spc, segment_cutoff=1):
    log.notice('Saving significant peaks.')
    cnt = 0
    for chrom in sorted(spc):
        segments = spc[chrom]
        segments.sort(key=operator.itemgetter(1)) # sort by start index
        for segment in segments:
            if segment.score >= segment_cutoff:
                cnt += 1
                of.write('\t'.join((chrom, str(segment.start),
                                    str(segment.end),
                                    str(segment.score))))
                of.write('\n')
    log.notice('Done. Wrote %d peaks.' % cnt)

def parse_chrom_filter(xs, prefix=''):
    if xs is None:
        return None
    if xs[0] == 'hg19':
        return hg19_chromfilter
    chroms = []
    for x in xs:
        if '-' in x:
            start, end = map(int, x.split('-'))
            chroms.extend(range(start, end + 1))
        else:
            chroms.append(x)
    return set('%s%s' % (prefix, x) for x in chroms)

def as_chrom_sizes(scores_per_chrom):
    chrom_sizes = {}
    npos_tot, ntot_tot = 0., 0.
    for chrom, spc in scores_per_chrom.items():
        ntot = len(spc)
        npos = sum(x.score for x in spc if x.score > 0)
        chrom_sizes[chrom] = (npos, ntot)
        npos_tot += npos
        ntot_tot += ntot
    r = float(npos_tot) / ntot_tot
    log.notice('\tPositive bin ratio is: %.2f' % r)
    assert r < 0.5, "Positive bin ratio must be less than 0.5"
    if r < .30:
        log.warn('Positive bin ratio is very low (%.2f)' % r)
    return chrom_sizes

def load_score_file(fname):
    '''ipython debug'''
    d = defaultdict(list)
    for line in open(fname):
        parts = line.split()
        x = bed(parts[0], int(parts[1]), int(parts[2]), float(parts[3]))
        d[x.chrom].append(x)
    return d
