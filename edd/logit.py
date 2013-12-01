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

def ci_for_df(odf, ci_min=0.25, min_reads_per_type=3, neg_score_scale=4):
    df = odf.copy()
    df['tot_reads'] = df.ip + df.input
    df['avg'] = df.ip / df.tot_reads.astype(float)
    __log_stats(df)
    z = 1.96
    p = df.avg.values.copy()
    const1 = p + (z**2)/(2*df.tot_reads)
    const2 = z * np.sqrt((p * (1-p) + z**2 / (4 * df.tot_reads)) / df.tot_reads)
    divisor = 1 +  z**2 / df.tot_reads
    df['ci_low'] = (const1 - const2) / divisor
    df['ci_high'] = (const1 + const2) / divisor
    df['ci_diff'] = df.ci_high - df.ci_low
    df['score'] = logit(df.ix[np.logical_and(df.tot_reads > 50, df.ci_diff < ci_min)].avg)
    median_neg, median_pos = get_medians(df.dropna())
    df.ix[np.logical_and(np.isnan(df.score),
                         df.avg > 0.5),
                         'score'] = median_pos
    df.ix[np.isnan(df.score), 'score'] = median_neg
    df.ix[df.score < 0, 'score'] *= neg_score_scale
    return df

def __log_stats(odf, at_least_lvl=48):
    df = odf.dropna()
    l = len(df.index) / 100.0 # gives percentage
    log.notice('cutoff to get at least %d: %.2f' % (at_least_lvl,
            np.percentile(df.avg, at_least_lvl)))
    negprct = (df.avg<= 0.5).sum() / l
    posprct = (df.avg> 0.5).sum() / l
    log.notice(('%.2f%% of bins are positive and %.2f%% of the bins are negative'
            % (posprct, negprct)))

