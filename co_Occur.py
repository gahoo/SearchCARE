#! /usr/bin/env python
#coding=utf-8
import sys
import os
import time
from zipfile import ZipFile, ZIP_DEFLATED
import cPickle
from math import log, exp, sqrt
from mpmath import loggamma

def co_Occurrence(Motif_LOC, NodeOutput='Node.txt', EdgeOutput='Edge.txt'):
    '''
    return co_Occur字典
    '''
    co_Occur={}
    Motif_stat={}
    NodeOutput=open(NodeOutput, 'w')
    EdgeOutput=open(EdgeOutput, 'w')
    #统计各motif的数量
    for motif in Motif_LOC.keys():
        Motif_stat[motif]=len(Motif_LOC[motif])
        NodeOutput.write(motif+'\t'+str(Motif_stat[motif])+'\n')
    #全集，拿来计算补集
    CompleteSet=set()
    for s in Motif_LOC.values():
        CompleteSet=CompleteSet|s
    LOC_count=len(CompleteSet)
    for motifa in Motif_LOC.keys():
        motifA=Motif_LOC.pop(motifa)
        motifA_=CompleteSet-motifA
        for motifb in Motif_LOC.keys():
            motifB=Motif_LOC[motifb]
            motifB_=CompleteSet-motifB
            keyname=motifa+'|'+motifb
            co_Occur[keyname]=[ \
                len(motifA&motifB), \
                len(motifA&motifB_), \
                len(motifA_&motifB), \
                len(motifA_&motifB_)]
            (a, b, c ,d)=[num for num in co_Occur[keyname]]
            (m, n, r, s)=[a+b, c+d, a+c, b+d]
            co_Occur[keyname].extend([m,n,r,s])
            #print keyname,(a, b, c ,d),(m, n, r, s)
            try:
                cooperationindex=(a*d-b*c)/sqrt(m*n*r*s)
                chi=float(LOC_count*pow(abs(a*d-b*c)-LOC_count/2,2))/(m*n*r*s)
                hypergeocdf=hypergeo_cdf(a,Motif_stat[motifa],Motif_stat[motifb],LOC_count)
                co_Occur[keyname].append(cooperationindex)
                co_Occur[keyname].append(chi)
                co_Occur[keyname].append(hypergeocdf)
            except ZeroDivisionError:
                print keyname,(a, b, c ,d),(m, n, r, s)
            EdgeOutput.write(motifa+'\t'+motifb+'\t'+'\t'.join([str(num) for num in co_Occur[keyname]])+'\n')
    return co_Occur

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

    assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
    assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
    assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
    assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)
    assert N-m >= n-X, 'There are more failures %i than unmarked items %i' % (N-m, n-X)

    s = 0
    for i in range(X, min(m,n)+1):
        s += max(gauss_hypergeom(i, n, m, N), 0.0)
    return min(max(s,0.0), 1)



if __name__ == "__main__":
    try:
        sys.argv[1]
    except:
        print "Tell me the filename please."
        sys.exit()
    if not os.path.exists(sys.argv[1]):
        raise Exception,"There is no file name %s" % sys.argv[1]
    basename=sys.argv[1].split('.')[0]
    nodefile="%s.Node" % basename
    edgefile="%s.Edge" % basename
    timestamp=time.time()
    #Motif_LOC = cPickle.load(open("%s.pkl" % basename,"rb"))
    zf = ZipFile('%s.zpkl' % basename, 'r')
    Motif_LOC = cPickle.loads(zf.open('Motif_LOC.pkl').read())
    zf.close()
    co_Occurrence(Motif_LOC,nodefile,edgefile)
    print time.time()-timestamp
