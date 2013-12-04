import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pa
from logbook import Logger

log = Logger(__name__)

def normal_opt_func(df):
    npos = df.npositive.values
    infoscore = df.information_score.values
    return npos * infoscore

class ScoreCutoff(object):
    """
    """

    def __init__(self, scores):
        """

        Arguments:
        - `scores`:
        """
        self.scores = scores
        self._ratio = None
        self.cutpoints = None
        self.cutpoint_scores = None
        self.df = None

    def __get_n_scores_uniformly(self, n):
        '''
        gets n (or close to n) values uniformly from all the
        score values.

        used to narrow search for cutoff value.
        '''
        ascores = np.unique(np.concatenate(self.scores))
        prcts = np.linspace(0, 1, n, endpoint=False)
        idx = np.unique((prcts * len(ascores)).astype(int))
        return ascores[idx]

    @classmethod
    def __information_score(cls, scores, pos_cutoff):
        '''
        Finds observed number of adjacent positive bins
        also: expected and log(observed/expected)

        bins needs to be binary
        todo -> should also work for score bins?
        '''
        nbins = sum(len(s) for s in scores)
        # there are N-1 adjacent bin pairs in an interval of N bins
        npairs = nbins - len(scores)
        npositive = sum((s > pos_cutoff).sum() for s in scores)
        r = float(npositive) / nbins
        expected = r**2 * npairs
        npositive, npositive_pairs = cls.__get_num_positive_bins_and_pairs(
                [x > pos_cutoff for x in scores])
        information_score = math.log((npositive_pairs + 1) / (expected + 1))
        return dict(score=(information_score * r),
                information_score=information_score, ratio=r,
                npositive=npositive, nbins=nbins, expected=expected,
                npositive_pairs=npositive_pairs, cutoff=pos_cutoff)

    def optimize(self, optf=normal_opt_func):
        '''
        Finds the cutpoint that gives the highest score.
        A cutpoints score is decided by the scoring function 
        `__information_score`.
        '''
        self.cutpoints = self.__get_n_scores_uniformly(300)

        log.notice('''\
        searching for optimal pos bin ratio lim between \
        %.3f and %.3f.''' % (
            self.cutpoints[0], self.cutpoints[-1]))

        self.df = pa.DataFrame([
            self.__information_score(self.scores, pos_cutoff)
            for pos_cutoff in self.cutpoints])
        xs = optf(self.df)
        xmin = xs[~np.isnan(xs)].min()
        xs[np.isnan(xs)] = xmin
        self.cutpoint_scores = xs
        self.max_idx = self.cutpoint_scores.argmax()
        self.lim_value = self.cutpoints[self.max_idx]
        return self

    def save_json(self, filename):
        self.df.to_csv(filename)


    @classmethod
    def from_chrom_scores(self, chrom_scores):
        scores = [np.array([x.score for x in xs])
                for xs in chrom_scores.values()]
        return ScoreCutoff(scores)

    def get_limit_score(self, ratio):
        ''' Returns the cutoff that yields the
        desired pos/neg ratio.
        '''
        bs = np.concatenate(self.scores)
        bs.sort()
        bs = bs[::-1]
        lim_score_idx = int(len(bs) * ratio)
        return bs[lim_score_idx]

##
## UTIL METHODS
##
    @classmethod
    def __get_num_positive_bins_and_pairs(cls, bins_by_chrom):
        npositive = npositive_pairs = 0
        # r = float(npos) / len(bins)
        # information_content = math.log((observed + 1) / float(expected + 1))
        for bins in bins_by_chrom:
            npositive_pairs += np.logical_and(bins[:-1], bins[1:]).sum()
            npositive += bins.sum()
        return npositive, npositive_pairs
    def __get_arrow_coords(self):
        ymax = self.cutpoint_scores[self.max_idx]
        x = self.cutpoints[self.max_idx]
        txt = 'x=%.2f, r=%.2f' % (x, self.ratio)
        return dict(s=txt, xy=(x, ymax),
                xytext=(self.cutpoints[len(self.cutpoints)/10.], ymax / 1.5),
                arrowprops=dict(facecolor='black', shrink=0.05))
    @property
    def ratio(self):
        if self._ratio is None:
            all_scores = np.concatenate(self.scores)
            orig_pos_ratio = (all_scores > 0).sum() / float(len(all_scores))
            self._ratio = (all_scores > self.lim_value).sum() / float(len(all_scores))
            log.notice('Original positive bin ratio is %.2f' % orig_pos_ratio)
            log.notice('Adjusted positive bin ratio is %.2f' % self._ratio)
        return self._ratio

    def plot_results(self, dst, title=None, annotate=False):
        plt.clf()
        plt.plot(self.cutpoints, self.cutpoint_scores)
        if title is not None:
            plt.title(title)
        if annotate:
            plt.annotate(**self.__get_arrow_coords())
        if dst is None:
            plt.show()
        else:
            plt.savefig(dst)
'''
 #  bin size optimizer
 def optimize(self, ic_power=1.0):
     # old
     xs = (np.sqrt(self.df.nbins.values) *
             self.df.ratio.values *
             self.df.information_score.values ** ic_power)

     # std of uniform (0, 1) dist
     ustd = 0.28865

     # skip 0 bins
     mccord_bin_std = 0.157352
     ad01_bin_std = 0.05811

     # as prct of ustd:
     mccord_spread = 0.545132618755
     ad01_spread = 0.201332132935

     xs = ((self.df.nbins.values ** (1. / (5 * (1 - spread)))) *
             self.df.ratio.values *
             self.df.information_score.values ** ic_power)
             '''
