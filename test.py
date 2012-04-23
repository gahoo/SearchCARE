from math import log, exp, sqrt
from mpmath import loggamma
from scipy import stats

def logchoose(ni, ki):
    try:
        lgn1 = loggamma(ni+1)
        lgk1 = loggamma(ki+1)
        lgnk1 = loggamma(ni-ki+1)
    except ValueError:
        #print ni,ki
        raise ValueError
    return lgn1 - (lgnk1 + lgk1)

def gauss_hypergeom(X, n, m, N):
    """Returns the probability of drawing X successes of m marked items
     in n draws from a bin of N total items."""

    assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
    assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
    assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
    assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)

    r1 = logchoose(m, X)
    try:
        r2 = logchoose(N-m, n-X)
    except ValueError:
        return 0
    r3 = logchoose(N,n)

    return exp(r1 + r2 - r3)

def hypergeo_cdf(X, n, m, N):
    '''
    '''

    assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
    assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
    assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
    assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)
    assert N-m >= n-X, 'There are more failures %i than unmarked items %i' % (N-m, n-X)

    s = 0
    for i in range(X, min(m,n)+1):
        s += max(gauss_hypergeom(i, n, m, N), 0.0)
    return min(max(s,0.0), 1)

def hypergeo_cdf_enrich(X, n, m, N):
    '''
    '''
    assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
    assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
    assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
    assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)
    assert N-m >= n-X, 'There are more failures %i than unmarked items %i' % (N-m, n-X)

    s = 0
    for i in range(0, X+1):
        s += max(gauss_hypergeom(i, n, m, N), 0.0)
    return min(max(s,0.0), 1)


test=[
#[11,187,11,285],
[5,6072,285,42086],
#[17,41829,285,42086],
[4,15622,9,42086],
[8,27694,9,42086],
[7,41929,9,42086],
[8,42035,9,42086],
[1,33646,9,42086],
[2,12159,9,42086],
[1,40181,9,42086],
[3,14328,9,42086],
[5,14180,9,42086],
[4,37197,9,42086],
[6,19036,9,42086],
[4,34155,9,42086],
[2,32381,9,42086],
[9,41723,9,42086],
[1,28455,9,42086],
[1,9574,9,42086],
[1,3461,9,42086],
[1,6464,9,42086],
[5,20635,9,42086],
[2,40581,9,42086],
[1,21608,9,42086],
[1,35994,9,42086],
[3,26854,9,42086],
[1,32459,9,42086],
[3,21039,9,42086],
[3,8127,9,42086],
[2,9895,9,42086],
[1,11731,9,42086],
[5,20509,9,42086],
[7,33898,9,42086],
[1,38975,9,42086],
[2,10204,9,42086],
[3,12506,9,42086],
[4,40007,9,42086],
[1,29578,9,42086],
[2,36444,9,42086],
[1,8238,9,42086],
[6,17083,9,42086],
[5,31226,9,42086],
[1,7234,9,42086],
[1,15338,9,42086],
[1,33801,9,42086],
[2,9656,9,42086],
[3,39003,9,42086],
[4,15622,9,42086],
[1,22390,9,42086],
]
#print hypergeo_cdf_enrich(5,6072,285,42086)
#print hypergeo_cdf_enrich(17,41829,285,42086)
for x,m,list_size,db_size in test:
    print 1-stats.hypergeom.cdf(x, db_size, m, list_size),1-hypergeo_cdf_enrich(x,list_size,m,db_size)