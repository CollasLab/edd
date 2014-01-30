import pandas as pa
from collections import namedtuple as __ntuple

Bedgraph = __ntuple('BedGraph', 'chrom start end score')

def read_bedgraph_as_dataframe(fname):
    return pa.read_table(fname, names=['chrom','start','end','score'])

def dataframe_as_list(df):
    return [Bedgraph(x['chrom'], x['start'], x['end'], x['score'])
            for _, x in df.iterrows()]

def read_bed_score_file(fname):
    xs = []
    with open(fname) as f:
        for line in f:
            chrom, start, end, score = line.split()
            xs.append(Bedgraph(chrom,int(start),int(end),float(score)))
    return xs

def scale_bedgraph_list(xs, scalef):
    return [Bedgraph(x.chrom, x.start, x.end, scalef(x.score))
            for x in xs]
