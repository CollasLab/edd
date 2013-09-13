def information_score_helper(bins):
    '''
    returns a score that tries to say something
    about coverage to quality.
    (we want to cover many high quality bins)
    '''
    npos = bins.sum()
    r = float(npos) / len(bins)
    expected = r**2 * (len(bins) - 1)
    observed = np.logical_and(bins[:-1], bins[1:]).sum()
    information_content = math.log((observed + 1) / float(expected + 1))
    return information_content * npos

def information_score(bins_per_chrom, pos_cutoff):
    return sum(information_score_helper(x > pos_cutoff) for x in bins_per_chrom)

def optimize_score_cutoff(scores):
    '''
    we require the cutoff to be between 0 and mean pos score
    '''
    all_scores = np.concatenate(scores)
    lim = all_scores[all_scores > 0].mean()
    log.notice('searching for optimal pos bin ratio lim between 0 and %.3f.' % lim)
    xs = np.linspace(0, lim, 70)
    ys = np.array([information_score(scores, pos_cutoff)
                   for pos_cutoff in xs])
    if verbose_dir is not None:
        opath = os.path.join(verbose_dir, 'information_content.png')
        plt.clf()
        plt.plot(xs, ys)
        plt.title('log2(obs / expected) * npos')
        plt.savefig(opath)
    lim_value = xs[ys.argmax()]
    orig_pos_ratio = (all_scores > 0).sum() / float(len(all_scores))
    pos_ratio = (all_scores > lim_value).sum() / float(len(all_scores))
    log.notice('Original positive bin ratio is %.2f' % orig_pos_ratio)
    log.notice('Adjusted positive bin ratio is %.2f' % pos_ratio)
    return lim_value
