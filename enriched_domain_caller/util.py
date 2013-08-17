from collections import namedtuple, defaultdict
import math
import operator

hg19_chromfilter = set(['chrY', 'chrX', 'chr13', 'chr12', 'chr11', 'chr10', 'chr17', 'chr16', 'chr15', 'chr14', 'chr19', 'chr18', 'chr22', 'chr20', 'chr21', 'chr7', 'chr6', 'chr5', 'chr4', 'chr3', 'chr2', 'chr1', 'chr9', 'chr8'])
bed = namedtuple('BedGraph', 'chrom start end score')
bin = namedtuple('Bin', 'chrom start end ip_count input_count')

def log2_score(pos, neg):
    if pos == 0 or neg == 0:
        return -5 # default negative for missing values
    return math.log(float(pos) / neg)
        

def obs_results(scores_per_chrom):
    segments_per_chrom = {}
    for k,v in scores_per_chrom.items():
        xs = max_segments([x.score for x in v])
        segments_per_chrom[k] = [bed(k, v[x.from_idx].start,
                                     v[x.to_idx].end, x.score) for x in xs]
    return segments_per_chrom

def filter_smaller_than_1sd_from_mean(d):
    scores = np.array([x.score for x in itertools.chain.from_iterable(d.values())])
    lim = scores.mean() + scores.std()
    return {k:[x for x in v if x.score > lim]
            for k, v in d.items()}


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
    return bins

def get_input_scale_factor(bins):
    ip_sum = input_sum = 0
    for x in bins:
        ip_sum += x.ip_count
        input_sum += x.input_count
    return float(ip_sum) / input_sum

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

def get_limit_score(bs, max_pos_bin_ratio):
    bs = list(bs)
    bs.sort(reverse=True, key=operator.attrgetter('score'))
    lim_score_idx = int(len(bs) * max_pos_bin_ratio)
    lim_score = bs[lim_score_idx].score
    return lim_score

def get_bedgraph_list(bs, lim_score):
    bs = []
    for x in bins:
        val = 1 if x.score >= lim_score else -1
        b = bed(x.chrom, x.start, x.end, val)
        bs.append(b)
    return b

def read_scores(bedgraph, legal_chroms, pos_bin_ratio):
    bins = read_counts(bedgraph, legal_chroms)
    input_scale_factor = get_input_scale_factor(bins)
    nbins = normalize_bins(bins, input_scale_factor)
    sbins = score_bins(nbins, ci_lower_bound)
    lim_score = get_limit_score(sbins, pos_bin_ratio)

    chromd = defaultdict(list)
    for x in binary_bins:
        bin_score = 1 if x.score >= lim_score else -1
        chromd[x.chrom].append(bed(x.chrom, x.start, x.end, bin_score))

    return chromd

def write_segments(of, spc, segment_cutoff=1):
    print 'Saving significant peaks.',
    for chrom in sorted(spc):
        segments = spc[chrom]
        segments.sort(key=operator.itemgetter(1)) # sort by start index
        for segment in segments:
            if segment.score >= segment_cutoff:
                of.write('\t'.join((chrom, str(segment.start),
                                    str(segment.end),
                                    str(segment.score))))
                of.write('\n')
    print 'Done'

def parse_chrom_filter(xs, prefix=''):
    if xs is None:
        return None
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
        npos = sum(x.score for x in spc if x.score == 1)
        chrom_sizes[chrom] = (npos, ntot)
        npos_tot += npos
        ntot_tot += ntot
    r = float(npos_tot) / ntot_tot
    print '\tPositive bin ratio is: %.2f' % r
    assert r < 0.5, "Positive bin ratio must be less than 0.5"
    return chrom_sizes
