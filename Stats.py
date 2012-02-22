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
    '''
    care=CAREdb(dbname)
    (SeqName_id,Motifs)=load2RAM(care,SeqName,Motifs)
    getDist(care,Motifs)
    (Motifs_SeqName,Motif_stat)=countSeq(care,Motifs)
    coOccur(care,SeqName_id,Motifs,Motifs_SeqName,Motif_stat)

def load2RAM(care,SeqName,Motifs):
    '''
    将指定部分数据库载入内存备用
    '''
    if SeqName:
        SeqName_id=[care.SeqName[sn] for sn in SeqName]
    else:
        SeqName_id=None
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
    care.cache('Scanned',['REF_SeqName','REF_MotifSeq','start'],SeqName=SeqName_id,MotifSeq=MotifSeq_id)
    return (SeqName_id,Motifs)

def getDist(care,Motifs):
    '''
    获得各个Motif在启动子上的分布情况
    '''
    motif_dist={}
    for motif in Motifs:
        motif_id=care.Motif[motif]
        SQL="""
        SELECT start
        FROM cache.Scanned
        WHERE REF_MotifSeq IN (
                SELECT REF_MotifSeq
                FROM cache.Instance
                WHERE REF_Motif=%s)
        """ % motif_id
        care.cur.execute(SQL)
        motif_dist[motif]=[rs[0] for rs in care.cur.fetchall()]
    zf = ZipFile("%s.zpkl" % basename, 'w', ZIP_DEFLATED)
    zf.writestr('motif_dist.pkl', cPickle.dumps(motif_dist))
    zf.close()

def countSeq(care,Motifs):
    '''
    统计Motif对应的SeqName_id，及其数量
    '''
    Motifs_SeqName={}
    Motif_stat={}
    for motif in Motifs:
        motif_id=care.Motif[motif]
        SQL="""
        SELECT REF_SeqName
        FROM cache.Scanned
        WHERE REF_MotifSeq IN (
                SELECT REF_MotifSeq
                FROM cache.Instance
                WHERE REF_Motif=%s)
        """ % motif_id
        care.cur.execute(SQL)
        Motifs_SeqName[motif_id]=set([rs[0] for rs in care.cur.fetchall()])
        Motif_stat[motif_id]=len(Motifs_SeqName[motif_id])
        Motifs_SeqName[motif_id]
        print motif,motif_id,Motif_stat[motif_id]
    return (Motifs_SeqName,Motif_stat)

def distance():
    return

def coOccur(care,SeqName_id,Motifs,Motifs_SeqName,Motif_stat,p=None):
    '''
    共现的Motif
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
    if SeqName_id:
        Seq_num=len(SeqName_id)
    else:
        Seq_num=len(care.SeqName)
    Motifs=dict([(k,None) for k in Motifs])
    m_total_num=len(Motifs)
    if not p:
        p=0.05
        #p=1.0/(m_total_num*(m_total_num-1))
    for motifa in Motifs.keys():
        motifA=Motifs_SeqName[care.Motif[motifa]]
        node.append("%s\t%s" % (motifa,Motif_stat[care.Motif[motifa]]))
        print len(Motifs)
        Motifs.pop(motifa)
        for motifb in Motifs:
            motifB=Motifs_SeqName[care.Motif[motifb]]
            m_ab="%s|%s" % (motifa,motifb)
            co_Occur[m_ab]=motifA&motifB
            co_num=len(motifA&motifB)
            a_num=Motif_stat[care.Motif[motifa]]
            b_num=Motif_stat[care.Motif[motifb]]
            hypergeocdf=hypergeo_cdf(co_num,a_num,b_num,Seq_num)
            if hypergeocdf<=p:
                row=[motifa,motifb,co_num,hypergeocdf]
                row="\t".join([str(col) for col in row])
                #print row
                edge.append(row)
    NodeOutput.write("\n".join(node))
    EdgeOutput.write("\n".join(edge))

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
    basename=sys.argv[1].split('.')[0]
    #nodefile="%s.Node" % basename
    #edgefile="%s.Edge" % basename
    timestamp=time.time()
    #Motif_LOC = cPickle.load(open("%s.pkl" % basename,"rb"))
    #zf = ZipFile('%s.zpkl' % basename, 'r')
    #Motif_LOC = cPickle.loads(zf.open('Motif_LOC.pkl').read())
    #zf.close()
    #co_Occurrence(Motif_LOC,nodefile,edgefile)
    ['Os12t0182500-01','Os12t0186600-03','Os12t0186901-00','Os12t0197150-00','Os12t0205633-00']
    ['HSE','Chs-CMA3','AAGAA-motif','CATTAT-motif','CTAG-motif']
    loc=loadLOC(sys.argv[2])
    Stats(sys.argv[1], \
    None, \
    None,)
    print time.time()-timestamp
