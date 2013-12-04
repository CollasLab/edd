import numpy as np
from logbook import Logger
log = Logger(__name__)

def logit(xs):
    return np.log(xs) - np.log(1 - xs)

def get_medians(df):
    neg = np.median(df.ix[df.score < 0].score)
    pos = np.median(df.ix[df.score > 0].score)
    #print df
    #print neg, pos
    return neg, pos

def get_ci_intervals(p, tot_reads):
    z = 1.96
    const1 = p + (z**2)/(2*tot_reads)
    const2 = z * np.sqrt((p * (1-p) + z**2 / (4 * tot_reads)) / tot_reads)
    divisor = 1 +  z**2 / tot_reads
    ci_low = (const1 - const2) / divisor
    ci_high = (const1 + const2) / divisor
    return ci_low, ci_high


def ci_for_df(odf, ci_min=0.25, min_reads_per_type=3, neg_score_scale=4, max_pos_prct=None):
    df = odf.copy()
    df['tot_reads'] = df.ip + df.input
    df['avg'] = df.ip / df.tot_reads.astype(float)
    df['ci_low'], df['ci_high'] = get_ci_intervals(df.avg, df.tot_reads)
    df['ci_diff'] = df.ci_high - df.ci_low
    df['score'] = logit(df.ix[np.logical_and(df.tot_reads > 50, 
        df.ci_diff < ci_min)].avg)
    median_neg, median_pos = get_medians(df.dropna())
    # positive scores led to some problems, when in doubt, be conservative
    #df.ix[np.logical_and(np.isnan(df.score),
    #                     df.avg > 0.5),
    #                     'score'] = median_pos
    df.ix[np.isnan(df.score), 'score'] = median_neg
    df.ix[df.score < 0, 'score'] *= neg_score_scale
    return df
