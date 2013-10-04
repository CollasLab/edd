import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pa
from logbook import Logger

log = Logger(__name__)



class ScoreCutoff(object):
    """
    """

    def __init__(self, scores, pos_bins_as_fraction=False):
        """

        Arguments:
        - `scores`:
        """
        self.scores = scores
        self.min_score = min(x.min() for x in self.scores)
        self.max_score = max(x.max() for x in self.scores)
        self._ratio = None
        self.pos_bins_as_fraction = pos_bins_as_fraction
        self.xs = None
        self.ys = None

    def __get_n_scores_uniformly(self, n):
        '''
        gets n (or close to n) values uniformly from all the
        score values.
        '''
        ascores = np.unique(np.concatenate(self.scores))
        prcts = np.linspace(0, 1, n, endpoint=False)
        idx = np.unique((prcts * len(ascores)).astype(int))
        return ascores[idx]

    def optimize(self, ic_power=1.0):
        log.notice('searching for optimal pos bin ratio lim between %.3f and %.3f.' % (
            self.min_score, self.max_score))

        self.xs = self.__get_n_scores_uniformly(300)
        self.df = pa.DataFrame([self.__information_score(self.scores, pos_cutoff)
                                for pos_cutoff in self.xs])
        xs = (np.sqrt(self.df.nbins.values) *
              self.df.ratio.values *
              self.df.information_score.values ** ic_power)
        xmin = xs[~np.isnan(xs)].min()
        # import ipdb
        # ipdb.set_trace()
        xs[np.isnan(xs)] = xmin
        self.ys =xs
        self.max_idx = self.ys.argmax()
        self.lim_value = self.xs[self.max_idx]
        return self

    def save_json(self, filename):
        self.df.to_csv(filename)


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
        plt.plot(self.xs, self.ys)
        if title is not None:
            plt.title(title)
        if annotate:
            plt.annotate(**self.__get_arrow_coords())
        plt.savefig(dst)

    def __get_arrow_coords(self):
        ymax = self.ys[self.max_idx]
        x = self.xs[self.max_idx]
        txt = 'x=%.2f, r=%.2f' % (x, self.ratio)
        return dict(s=txt, xy=(x, ymax),
                     xytext=(self.xs[len(self.xs)/10.], ymax / 1.5),
                     arrowprops=dict(facecolor='black', shrink=0.05))

    @classmethod
    def from_chrom_scores(self, chrom_scores, pos_bins_as_count=False):
        scores = [np.array([x.score for x in xs])
                  for xs in chrom_scores.values()]
        return ScoreCutoff(scores, pos_bins_as_fraction=(not pos_bins_as_count))

    def get_limit_score(self, ratio):
        bs = np.concatenate(self.scores)
        bs.sort()
        bs = bs[::-1]
        lim_score_idx = int(len(bs) * ratio)
        return bs[lim_score_idx]

    @classmethod
    def __get_num_positive_bins_and_pairs(cls, bins_by_chrom):
        npositive = npositive_pairs = 0
        # r = float(npos) / len(bins)
        # information_content = math.log((observed + 1) / float(expected + 1))
        for bins in bins_by_chrom:
            npositive_pairs += np.logical_and(bins[:-1], bins[1:]).sum()
            npositive += bins.sum()
        return npositive, npositive_pairs

    @classmethod
    def __information_score(cls, scores, pos_cutoff):
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
