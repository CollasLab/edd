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


def ci_for_df(odf, ci_min=0.25, pscore_lim=10, neg_score_scale=4,
        extrapolate_low_info_bins=True):
    df = odf.copy()
    df['tot_reads'] = df.ip + df.input
    df['avg'] = df.ip / df.tot_reads.astype(float)
    df['ci_low'], df['ci_high'] = get_ci_intervals(df.avg, df.tot_reads)
    df['ci_diff'] = df.ci_high - df.ci_low
    # we assume that the sample mean is normally distributed
    # if equation below is > 10. If so, we can compute the
    # 95% confidence interval
    norm_sample_mean = np.minimum(df.avg, 1 - df.avg) * df.tot_reads > pscore_lim
    small_CI = df.ci_diff < ci_min
    scorable_bins = np.logical_and(norm_sample_mean, small_CI)
    df['score'] = logit(df.ix[scorable_bins].avg)
    df.ix[df.score < 0, 'score'] *= neg_score_scale
    if extrapolate_low_info_bins:
        median_neg, median_pos = get_medians(df.dropna())
        # positive scores led to some problems, when in doubt, be conservative
        #df.ix[np.logical_and(np.isnan(df.score),
        #                     df.avg > 0.5),
        #                     'score'] = median_pos
        df.ix[np.isnan(df.score), 'score'] = median_neg
    return df

def get_nib_ratio(df, max_ci_diff):
    df = df.copy()
    df['tot'] = df.ip + df.input
    df['avg'] = df.ip / df.tot.astype(float)
    nbins_with_reads = (df.tot > 0).sum()
    req = np.minimum(df.avg, 1 - df.avg) * df.tot > 10
    df = df.ix[req]
    ci_low, ci_high = get_ci_intervals(df.avg, df.tot)
    nbins_ok = (ci_high - ci_low < max_ci_diff).sum()
    return float(nbins_with_reads - nbins_ok) / nbins_with_reads
