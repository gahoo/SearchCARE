#! /usr/bin/env python
#coding=utf-8
import sys
import os
import time
import cPickle
from math import log, exp, sqrt
from mpmath import loggamma
from zipfile import ZipFile, ZIP_DEFLATED
from CAREdb import CAREdb

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

def Stats(dbname, SeqName=None, Motifs=None):
    '''
    统计Motif的数目，分布，共现情况
    motif_dist[motif][MotifSeq]={'+':dist,'-':dist}
    Motif_SeqName[motif][MotifSeq]={'+':SeqName,'-':SeqName}
    '''
    care=CAREdb(dbname)
    basename=care.dbfile.split('.')[0]
    (REF_SeqName,Motifs)=loadDB2RAM(care,SeqName,Motifs)
    id2MotifSeq=dict([(v,k) for k,v in care.MotifSeq.items()])
    motif_dist={}
    Motif_SeqName={}
    for motif in Motifs:
        REF_MotifSeq=getREF_MotifSeq(care,motif)
        motif_dist[motif]={}
        Motif_SeqName[motif]={}
        for MotifSeq_id in REF_MotifSeq:
            MotifSeq=id2MotifSeq[MotifSeq_id]
            motif_dist[motif][MotifSeq]=getDist(care,MotifSeq_id)
            Motif_SeqName[motif][MotifSeq]=getREF_SeqName(care,MotifSeq_id)
    coOccur(care,REF_SeqName,Motifs,Motif_SeqName)
    zf = ZipFile("%s.zpkl" % basename, 'w', ZIP_DEFLATED)
    zf.writestr('motif_dist.pkl', cPickle.dumps(motif_dist))
    zf.writestr('Motif_SeqName.pkl', cPickle.dumps(Motif_SeqName))
    zf.close()

def loadDB2RAM(care,SeqName,Motifs):
    '''
    将指定部分数据库载入内存备用
    '''
    if SeqName:
        REF_SeqName=[care.SeqName[sn] for sn in SeqName]
    else:
        REF_SeqName=None
    if Motifs:
        MotifSeq_id=[]
        for motif in Motifs:
            SQL="SELECT REF_MotifSeq FROM Instance WHERE REF_Motif=%s" % care.Motif[motif]
            care.cur.execute(SQL)
            MotifSeq_id.extend([rs[0] for rs in care.cur.fetchall()])
        MotifSeq_id=list(set(MotifSeq_id))
    else:
        MotifSeq_id=None
        Motifs=care.Motif.keys()
    print "Load to Cache"
    care.cache('Instance')
    care.cache('Scanned',['REF_SeqName','REF_MotifSeq','start','stop','strand'],SeqName=REF_SeqName,MotifSeq=MotifSeq_id)
    return (REF_SeqName,Motifs)

def getREF_MotifSeq(care,motif):
    motif_id=care.Motif[motif]
    SQL="""
    SELECT DISTINCT REF_MotifSeq
    FROM cache.Instance
    WHERE REF_Motif=%s
    """ % motif_id
    care.cur.execute(SQL)
    REF_MotifSeq=[rs[0] for rs in care.cur.fetchall()]
    return REF_MotifSeq

def getDist(care,MotifSeq_id,strand='both'):
    '''
    获得MotifSeq在启动子上的分布情况
    '''
    if strand=='both':
        WHERE=""
    elif strand=='+':
        WHERE="AND strand=1"
    elif strand=='-':
        WHERE="AND strand=0"
    SQL="""
    SELECT start,stop,strand
    FROM cache.Scanned
    WHERE REF_MotifSeq=%s
        %s
    """ % (MotifSeq_id,WHERE)
    care.cur.execute(SQL)
    dist={'+':[],'-':[]}
    for start,stop,strand in care.cur.fetchall():
        if strand==1:
            dist['+'].append(start)
        else:
            dist['-'].append(stop)
    return dist

def getREF_SeqName(care,MotifSeq_id,strand='both'):
    '''
    统计REF_MotifSeq对应的REF_SeqName，及其数量
    TODO 细化到具体的序列
    '''
    if strand=='both':
        WHERE=""
    elif strand=='+':
        WHERE="AND strand=1"
    elif strand=='-':
        WHERE="AND strand=0"
    SQL="""
    SELECT REF_SeqName,strand
    FROM cache.Scanned
    WHERE REF_MotifSeq=%s
        %s
    """ % (MotifSeq_id,WHERE)
    care.cur.execute(SQL)
    SeqName={'+':[],'-':[]}
    for REF_SeqName,strand in care.cur.fetchall():
        if strand==1:
            SeqName['+'].append(REF_SeqName)
        else:
            SeqName['-'].append(REF_SeqName)
    return SeqName

def coOccur(care,REF_SeqName,Motifs,Motif_SeqName,p=None):
    '''
    共现的Motif
    TODO 细化到具体的序列
    TODO 考虑距离问题，若是距离大于某个值，就不算?
    '''
    basename=care.dbfile.split('.')[0]
    nodefile="%s.Node" % basename
    edgefile="%s.Edge" % basename
    NodeOutput=open(nodefile, 'w')
    EdgeOutput=open(edgefile, 'w')
    node=[]
    edge=[]
    co_Occur={}
    print "-------"
    if REF_SeqName:
        Seq_num=len(REF_SeqName)
    else:
        Seq_num=len(care.SeqName)
    Motifs=dict([(k,None) for k in Motifs])
    m_total_num=len(Motifs)
    if not p:
        p=0.05
        #p=1.0/(m_total_num*(m_total_num-1))
    for motifa in Motifs.keys():
        motifA=mergeSeqName(Motif_SeqName[motifa])
        node.append("%s\t%s" % (motifa,len(motifA)))
        print len(Motifs)
        Motifs.pop(motifa)
        for motifb in Motifs:
            motifB=mergeSeqName(Motif_SeqName[motifb])
            m_ab="%s|%s" % (motifa,motifb)
            co_Occur[m_ab]=motifA&motifB
            co_num=len(motifA&motifB)
            #此处会重复统计SeqName的个数，最好能够用一个变量直接记起来，直接查询
            a_num=len(motifA)
            b_num=len(motifB)
            hypergeocdf=hypergeo_cdf(co_num,a_num,b_num,Seq_num)
            if hypergeocdf<=p:
                row=[motifa,motifb,co_num,hypergeocdf]
                row="\t".join([str(col) for col in row])
                #print row
                edge.append(row)
    NodeOutput.write("\n".join(node))
    EdgeOutput.write("\n".join(edge))

def mergeSeqName(motif_dict):
    '''
    motif_dict={MotifSeq1:{'+':[],'-':[]},}
    '''
    Seq_id=[]
    [Seq_id.extend(v['+']+v['-']) for v in motif_dict.values()]
    return set(Seq_id)

def distance():
    '''
    统计两个motif间的距离
    TODO 细化到具体的序列
    '''
    return


def loadLOC(filename):
    return [line.strip() for line in open(filename).readlines()]

if __name__ == "__main__":
    try:
        sys.argv[1]
    except:
        print "Tell me the filename please."
        sys.exit()
    if not os.path.exists(sys.argv[1]):
        raise Exception,"There is no file name %s" % sys.argv[1]
    timestamp=time.time()
    ['Os12t0182500-01','Os12t0186600-03','Os12t0186901-00','Os12t0197150-00','Os12t0205633-00']
    ['HSE','Chs-CMA3','AAGAA-motif','CATTAT-motif','CTAG-motif']
    loc=loadLOC(sys.argv[2])
    Stats(sys.argv[1], \
    None, \
    None,)
    print time.time()-timestamp
