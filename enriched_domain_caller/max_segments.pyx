from collections import namedtuple
from libc.stdlib cimport malloc, free, realloc
cimport cython

max_segment_struct = namedtuple('MaxSegment', 'score from_idx to_idx')

cdef struct segment:
    int I
    int L
    int R
    int lidx
    int ridx

cdef int consider(segment *s, segment *buf, int k):
    cdef int j = k - 1
    while (j >= 0 and buf[j].L >= s.L):
        j -= 1
    if (j == -1 or buf[j].R >= s.R):
        buf[k] = s[0]
        return k + 1
    else:
        s.I = s.R - buf[j].L
        s.L = buf[j].L
        s.lidx = buf[j].lidx
        s.ridx = s.ridx
        return consider(s, buf, j)

cdef enum:
    N = 1000

cdef int max_segments_impl(int *xs, int len_xs,
                           segment **buf, int len_buf):
    cdef int csum = 0
    cdef int k = 0
    cdef int i
    cdef segment s
    cdef int x
    cdef int Lk
    cdef int Rk
    for i in range(len_xs):
        x = xs[i]
        Lk = csum
        Rk = csum + x
        csum += x;
        if (k == len_buf):
            len_buf *= 2
            if (len_buf < N):
                # if len_xs is initially zero
              len_buf = N
            buf[0] = <segment *> realloc(buf[0], len_buf * sizeof(segment))
        if (x > 0):
            s = segment(x, Lk, Rk, i, i)
            k = consider(&s, buf[0], k)
    return k

def max_segments(xs):
    cdef int *cxs = <int*> malloc(len(xs)*cython.sizeof(int))
    cdef segment *buf = NULL
    for i, x in enumerate(xs):
        cxs[i] = x
    k = max_segments_impl(cxs, len(xs), &buf, 0)
    free(cxs)
    res = []
    for i in range(k):
        y = max_segment_struct(buf[i].I, buf[i].lidx, buf[i].ridx)
        res.append(y)
    free(buf)
    return res
