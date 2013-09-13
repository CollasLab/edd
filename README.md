edc
===

Enriched Domain Caller


TODO
====

* Monte Carlo simulation must retain bin scores (now it assumes that all bins are either +1 or -1).
* Issue warning if bin size is too small / quality is too low (information_score * nposbins graph should be concave)
* find good reasoning behind lower estimated FDR cutoff
* is there some background noise on chrY for females? should it be excluded if read count is below a certain threshold?
* DONE - mask large gaps
* Add tool for optimize cutoff plot (concave stuff), so that different bin sizes can be checked
* Require that input data is normalized (perhaps with a flag that does simple normalization)