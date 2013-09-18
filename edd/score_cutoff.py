import math

import numpy as np
import matplotlib.pyplot as plt
from logbook import Logger

log = Logger(__name__)



class ScoreCutoff(object):
    """
    """

    def __init__(self, scores):
        """

        Arguments:
        - `scores`:
        """
        self.scores = scores
        self.min_score = min(x.min() for x in self.scores)
        self.max_score = max(x.max() for x in self.scores)
        self._ratio = None

    def optimize(self):
        log.notice('searching for optimal pos bin ratio lim between %.3f and %.3f.' % (
            self.min_score, self.max_score))
        self.xs = np.linspace(self.min_score, self.max_score, 1000)
        self.ys = self.__information_score_for_range()
        self.max_idx = self.ys.argmax()
        self.lim_value = self.xs[self.max_idx]
        return self

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
    def from_chrom_scores(self, chrom_scores):
        scores = [np.array([x.score for x in xs])
                  for xs in chrom_scores.values()]
        return ScoreCutoff(scores)

    def get_limit_score(self, ratio):
        bs = np.concatenate(self.scores)
        bs.sort()
        bs = bs[::-1]
        lim_score_idx = int(len(bs) * ratio)
        return bs[lim_score_idx]

    @classmethod
    def information_score_helper(cls, bins):
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

    def __information_score(self, pos_cutoff):
        return sum(self.information_score_helper(x > pos_cutoff)
                   for x in self.scores)

    def __information_score_for_range(self):
        return np.array([self.__information_score(pos_cutoff)
                         for pos_cutoff in self.xs])
