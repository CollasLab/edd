import collections
import operator

import logbook
import pybedtools


log = logbook.Logger(__name__)

class Gap(object):

    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end

    def bigger_than(self, x):
        '''
        returns true if gap comes after param x (no overlap)
        '''
        return x.end <= self.start

    def overlaps(self, x):
        '''
        returns true if the two intervals overlaps
        '''
        return (self.start <= x.start < self.end or
                self.start < x.end <= self.end)

    def __repr__(self):
        return 'Gap(%s, %d, %d)' % (self.chrom, self.start, self.end)

def read_gap_file(path, drop_smaller_than):
    res = []
    src = pybedtools.BedTool(path)
    filtered_src = src.merge().filter(lambda x: (x.end - x.start) >= drop_smaller_than)
    log.notice('Gap file read. Total coverage: %.2fMB' % (
        src.total_coverage() / 1e6))
    for e in filtered_src:
        x = Gap(e.chrom, e.start, e.end)
        res.append(x)
    tot = sum(x.end - x.start for x in res)
    log.notice('Removing gaps smaller than %.2fMB. Total coverage after filtering: %.2fMB' % (
        (drop_smaller_than / 1e6, tot / 1e6)))
    return res

def split_on_gaps(scores_per_chrom, gaps):
    dg = collections.defaultdict(list)
    revdict = {}
    d = {}
    for g in gaps:
        dg[g.chrom].append(g)

    for chrom, gaps in dg.items():
        if not chrom in scores_per_chrom:
            continue
        gaps.sort(key=operator.attrgetter('start'))
        groups = []
        cur_grp = []
        cur_gap = 0
        ndropped = 0
        bins = iter(scores_per_chrom[chrom])
        try:
            while True:
                x = bins.next()
                if gaps[cur_gap].bigger_than(x):
                    cur_grp.append(x)
                elif gaps[cur_gap].overlaps(x):
                    if len(cur_grp) > 0:
                        groups.append(cur_grp)
                        cur_grp = []
                    ndropped += 1
                else:
                    cur_gap += 1
                    cur_grp.append(x)
                    if len(gaps) == cur_gap:
                        cur_grp.extend(bins)
                        break
                    assert not gaps[cur_gap].overlaps(x)
        except StopIteration:
            pass
        if len(cur_grp) > 0:
            groups.append(cur_grp)
        names = ['%s_%d' % (chrom, i) for i in range(len(groups))]
        d.update(zip(names, groups))
        for name in names:
            revdict[name] = chrom
    for chrom in (scores_per_chrom.viewkeys() - dg.viewkeys()):
        d[chrom] = scores_per_chrom[chrom]
    return d, revdict

def join_gaps(res, gapped_chroms_to_chrom):
    d = collections.defaultdict(list)
    for k, xs in res.items():
        d[gapped_chroms_to_chrom.get(k, k)].extend(xs)
    for xs in d.values():
        xs.sort(key=operator.attrgetter('start'))
    return d

