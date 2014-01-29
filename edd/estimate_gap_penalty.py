import logit, util
import pandas as pa
import StringIO
import tempfile
import os

def count_stats(xs):
    '''xs is a bedtool instance where the name field holds the bin score'''
    stats = {'DIB': 0, 'EIB': 0}
    for x in xs:
        if float(x.name) > 0:
            stats['EIB'] += 1
        else:
            stats['DIB'] += 1
    return stats

def estimate_gap_penalty(odf, nprocs, mc_trials=100, outfile_path=None):
    binscore_df = logit.ci_for_df(odf, neg_score_scale=1)
    bedgraph_path = tempfile.mktemp()
    util.save_bin_score_file(binscore_df, bs_path)
    bg = BedTool(bedgraph_path)
    xs = []
    
    for neg_score_scale in range(2,20):
        df = edd.logit.ci_for_df(odf, neg_score_scale=neg_score_scale)
        gb = edd.df_as_bins(df, args.gap_file.name)
        max_bin_score = df.score.max()
        observed_result = gb.max_segments(filter_trivial=max_bin_score)
        mc_res = edd.algorithm.MonteCarlo.run_simulation(gb.chrom_scores, 
                                                         niter=mc_trials, nprocs=nprocs)
        tester = edd.algorithm.max_segments.IntervalTest(observed_result, mc_res)
        segments = [segment for (segment, pval) in tester.pvalues() if pval < 0.05]
        peaks_sb = StringIO.StringIO()
        tester.segments_to_bedstream(segments, peaks_sb)
        peaks = BedTool(peaks_sb.get_value(), from_string=True)
        d = count_stats(bg.intersect(peaks))
        d['gap-penalty'] = neg_score_scale
        xs.append(d)
    genome_wide_stats = count_stats(bg)
    df = pa.DataFrame(xs)
    df['peak_EIB_ratio'] = df.EIB / (df.EIB + df.DIB).astype(float)
    df['global_EIB_coverage'] = df.EIB / float(genome_wide_stats['EIB'])
    df['score'] = df.peak_EIB_ratio**5 * df.global_EIB_coverage
    df.sort('gap-penalty', inplace=True)
    if outfile_path:
        df.to_csv(outfile_path, index=False)
    # TODO perform extra gap_penalty tests for best score +- 0.5
    os.remove(bedgraph_path)
    return df.ix[df.score.argmax()]['gap-penalty']
