import util
from max_segments import max_segments
import unittest
from nose.tools import assert_almost_equal

def test_trivial():
    xs = max_segments([1,1,1])
    x, = xs
    assert x.score == 3

def test_trivial_gap():
    xs = max_segments([1,-1,1])
    x,y = xs
    assert x.score == 1
    assert y.score == 1

def test_merged_gap():
    xs = max_segments([1,1,-1,1,1])
    x, = xs
    assert x.score == 3

def test_largescore_trivial():
    xs = max_segments([3,2])
    x, = xs
    assert x.score == 5

def test_largescore_wnegative():
    xs = max_segments([3,-2, 4])
    x, = xs
    assert x.score == 5

def test_largescore_float():
    xs = max_segments([3.5,2.5])
    x, = xs
    assert x.score == 6

def test_largescore_negfloat():
    xs = max_segments([3.5,-0.5,2.5])
    x, = xs
    assert_almost_equal(x.score, 5.5)
