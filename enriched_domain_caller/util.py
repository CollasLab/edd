from collections import namedtuple, defaultdict
import operator

hg19_chromfilter = set(['chrY', 'chrX', 'chr13', 'chr12', 'chr11', 'chr10', 'chr17', 'chr16', 'chr15', 'chr14', 'chr19', 'chr18', 'chr22', 'chr20', 'chr21', 'chr7', 'chr6', 'chr5', 'chr4', 'chr3', 'chr2', 'chr1', 'chr9', 'chr8'])
bed = namedtuple('BedGraph', 'chrom start end score')

def read_counts(bedgraph, legal_chroms):
    bin = namedtuple('Bin', 'chrom start end ip_count input_count')
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

def get_limit_score(bins, max_pos_bin_ratio, input_scale_factor):
    bs = []
    for x in bins:
        if x.ip_count > 0 and x.input_count > 0:
            score = x.ip_count / (x.ip_count + x.input_count * input_scale_factor)
        else:
             score = 0
        bs.append(bed(x.chrom, x.start, x.end, score))
    bs.sort(reverse=True, key=operator.attrgetter('score'))
    lim_score_idx = int(len(bs) * max_pos_bin_ratio)
    while bs[lim_score_idx].score <= 0.5:
        lim_score_idx -= 1
        assert lim_score_idx >= 0
    lim_score = bs[lim_score_idx].score
    return lim_score, bs

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
    lim_score, bgs = get_limit_score(bins, pos_bin_ratio, input_scale_factor)
    chromd = defaultdict(list)
    for x in bgs:
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
