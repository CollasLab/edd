def count_stats(xs):
    stats = {'DIB': 0, 'EIB': 0}
    for x in xs:
        if float(x.name) > 0:
            stats['EIB'] += 1
        else:
            stats['DIB'] += 1
    return stats

def stats_for_peaks(peak_path, bedgraph_path):
    ' TODO del? '
    peaks = BedTool(peak_path)
    bg = BedTool(bedgraph_path)
    stats = {'name': peak_path,
             'neg-score': path_to_neg_score(os.path.basename(peak_path))}
    stats.update(count_stats(bg.intersect(peaks)))
    return stats
    
def df_stats(xs, genome_wide):
    df = pa.DataFrame(xs)
    del df['name']
    df['peak_EIB_ratio'] = df.EIB / (df.EIB + df.DIB).astype(float)
    df['global_EIB_coverage'] = df.EIB / float(genome_wide['EIB'])
    # not sure if the equation below is meaningful
    df['peak_EIB_weight_ratio'] = df.EIB_sum / (df.EIB_sum + df.DIB_sum).astype(float)
    df['global_EIB_weight_coverage'] = df.EIB_sum / float(genome_wide['EIB_sum'])
    
    return df.sort('neg-score').set_index('neg-score')

def read_df_stats(peak_paths, bg_path):
    return df_stats([stats_for_peaks(p,bg_path) 
                     for p in peak_paths],
                    count_stats())
    
def estimate_gap_penalty(odf, bedgraph_path, nprocs, mc_trials=100, outfile_path=None):
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
    if outfile_path is not None:
        df.to_csv(outfile_path, index=False)
    return df.ix[df.score.argmax()]['gap-penalty']
