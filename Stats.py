#! /usr/bin/env python
#coding=utf-8
import sys
import os
import time
import getopt
import cPickle
from scipy import stats
from sqlite3 import OperationalError
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
    '''
    用来计算共现的p
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
    用来计算富集的p
    N：数据库中基因总数
    m：数据库中特定motif数
    n：列表中基因总数
    x：列表中特定motif数
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

def Stats(dbname, SeqName=None, Motifs=None, p=None, output=None):
    '''
    统计Motif的数目，分布，共现情况
    motif_dist[motif][MotifSeq]={'+':dist,'-':dist}
    Motif_SeqName[motif][MotifSeq]={'+':SeqName,'-':SeqName}
    '''
    care=CAREdb(dbname)
    if not output:
        basename=care.dbfile.split('.')[0]
    else:
        basename=output
    checkMotif_Counts(care)
    (REF_SeqName,Motifs)=loadDB2RAM(care,SeqName,Motifs)
    id2MotifSeq=dict([(v,k) for k,v in care.MotifSeq.items()])
    motif_dist={}
    Motif_SeqName={}
    Motif_Enrich={}
    for motif in Motifs:
        REF_MotifSeq=getREF_MotifSeq(care,motif)
        motif_dist[motif]={}
        Motif_SeqName[motif]={}
        for MotifSeq_id in REF_MotifSeq:
            MotifSeq=id2MotifSeq[MotifSeq_id]
            motif_dist[motif][MotifSeq]=getDist(care,MotifSeq_id)
            Motif_SeqName[motif][MotifSeq]=getREF_SeqName(care,MotifSeq_id)
    MergedSeqName=getMerged(Motif_SeqName,'SeqName')
    Merged_dist=getMerged(motif_dist,'dist')
    SeqName_counts=countSeqName(MergedSeqName)
    zf = ZipFile("%s.zpkl" % basename, 'w', ZIP_DEFLATED)
    if SeqName:
        Enriched=enrichMent(care,SeqName_counts,list_size=len(REF_SeqName),output=output)
        zf.writestr('Enriched.pkl', cPickle.dumps(Enriched))
    zf.writestr('motif_dist.pkl', cPickle.dumps(motif_dist))
    zf.writestr('Merged_dist.pkl', cPickle.dumps(Merged_dist))
    zf.writestr('Motif_SeqName.pkl', cPickle.dumps(Motif_SeqName))
    zf.writestr('MergedSeqName.pkl', cPickle.dumps(MergedSeqName))
    zf.writestr('SeqName_counts.pkl', cPickle.dumps(SeqName_counts))
    co_Occur=coOccur(care,REF_SeqName,SeqName_counts,MergedSeqName,p,output)
    zf.writestr('co_Occur.pkl', cPickle.dumps(co_Occur))
    zf.close()

def loadDB2RAM(care,SeqName,Motifs):
    '''
    将指定部分数据库载入内存备用
    '''
    if SeqName:
        REF_SeqName=[care.SeqName[sn] for sn in SeqName if care.SeqName.has_key(sn)]
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

def enrichMent(care,SeqName_counts,list_size,output=None):
    '''
    计算Motif是否在输入列表中富集
    '''
    if not output:
        basename=care.dbfile.split('.')[0]
    else:
        basename=output
    db_size=len(care.SeqName)
    nodefile="%s.Node" % basename
    enrichNodefile="%s_enriched.Node" % basename
    Motif_Desc=care.loadDB2dict('Motif',1,2)
    node=[]
    enrichNode=[]
    Enriched={}
    for motif in SeqName_counts.keys():
        print motif
        REF_Motif=care.Motif[motif]
        x=SeqName_counts[motif]
        m=care.Motif_Counts[REF_Motif]
        if not x:
            print "not found"
            continue
        ER=(float(x)/list_size)/(float(m)/db_size)
        print x,m,list_size,db_size,ER,stats.hypergeom.cdf(x,db_size,m,list_size)
        try:
            p=1-hypergeo_cdf_enrich(x,list_size,m,db_size)
        except AssertionError,ERR:
            #17,285,41829,42086
            print ERR
            p=1
        print x,m,list_size,db_size,ER,p
        row=(motif,Motif_Desc[motif],x,m,list_size,db_size,ER,p)
        row=[str(r) for r in row]
        node.append("\t".join(row) )
        if p<0.05:
            Enriched[motif]={}
            Enriched[motif]['ER']=ER
            Enriched[motif]['p']=p
            enrichNode.append("\t".join(row) )
    NodeOutput=open(nodefile, 'w')
    NodeOutput.write("\n".join(node))
    enrichNodefile=open(enrichNodefile, 'w')
    enrichNodefile.write("\n".join(enrichNode))
    return Enriched

def coOccur(care,REF_SeqName,SeqName_counts,MergedSeqName,p=None,output=None):
    '''
    共现的Motif
    TODO 细化到具体的序列
    TODO 考虑距离问题，若是距离大于某个值，就不算?
    '''
    if not output:
        basename=care.dbfile.split('.')[0]
    else:
        basename=output
    edgefile="%s.Edge" % basename
    EdgeOutput=open(edgefile, 'w')
    edge=[]
    co_Occur={}
    print "-------"
    if REF_SeqName:
        Seq_num=len(REF_SeqName)
    else:
        Seq_num=len(care.SeqName)
    m_total_num=len(MergedSeqName)
    if not p:
        p=0.05
        #p=1.0/(m_total_num*(m_total_num-1))
    for motifa in MergedSeqName.keys():
        motifA=MergedSeqName[motifa]
        print len(MergedSeqName)
        MergedSeqName.pop(motifa)
        for motifb in MergedSeqName.keys():
            motifB=MergedSeqName[motifb]
            m_ab="%s|%s" % (motifa,motifb)
            co_Occur[m_ab]=motifA&motifB
            co_num=len(motifA&motifB)
            a_num=SeqName_counts[motifa]
            b_num=SeqName_counts[motifb]
            hypergeocdf=hypergeo_cdf(co_num,a_num,b_num,Seq_num)
            if hypergeocdf<=p:
                row=[motifa,motifb,co_num,a_num,b_num,hypergeocdf]
                row="\t".join([str(col) for col in row])
                #print row
                edge.append(row)
    EdgeOutput.write("\n".join(edge))
    return co_Occur

def getMerged(motif_dict,type):
    '''
    合并Motif_SeqName结果，即把：
    Motif_SeqName[motif][MotifSeq]={'+':SeqName,'-':SeqName}
    变成：
    Motif_SeqName[motif]=SeqName
    '''
    Merged={}
    if type=='dist':
        for motif in motif_dict.keys():
            Merged[motif]=merge(motif_dict[motif])
    elif type=='SeqName':
        for motif in motif_dict.keys():
            Merged[motif]=mergeSeqName(motif_dict[motif])
    return Merged

def merge(motif_dict):
    '''
    motif_dict={MotifSeq1:{'+':[],'-':[]},}
    '''
    Seq_id=[]
    [Seq_id.extend(v['+']+v['-']) for v in motif_dict.values()]
    return Seq_id

def mergeSeqName(motif_dict):
    '''
    motif_dict={MotifSeq1:{'+':[],'-':[]},}
    '''
    return set(merge(motif_dict))

def countSeqName(MergedSeqName):
    '''
    统计SeqName数目
    '''
    SeqName_counts={}
    for motif in MergedSeqName.keys():
        SeqName_counts[motif]=len(MergedSeqName[motif])
    return SeqName_counts

def distance():
    '''
    统计两个motif间的距离
    TODO 细化到具体的序列
    '''
    return

def loadList(filename):
    '''
    加载列表文件内容
    '''
    lst=list(set([line.strip() for line in open(filename).readlines()]))
    print "%s loaded" % len(lst)
    try:
        lst.remove('')
    except ValueError:
        pass
    return lst

def checkMotif_Counts(care):
    '''
    检查数据库内Motif_Counts是否有内容，没有就统计
    '''
    if not care.Motif_Counts:
        print "first stats"
        countAllSeqNameInDB(care)
        care.Motif_Counts=care.loadDB2dict('Motif_Counts',0,1)

def countAllSeqNameInDB(care):
    '''
    统计数据库内全部Motif对应的SeqName数目
    #TODO 必须注意，如果选取的序列长度不同应该重新进行计算，否则背景是不同的，会导致Enrichment的计算有误，可以考虑增加一列记录长度的字段
    '''
    SQL="""
    INSERT INTO Motif_Counts (REF_Motif,SeqName_counts)
    VALUES (%s,(
        SELECT COUNT(DISTINCT REF_SeqName)
        FROM Scanned
        WHERE REF_MotifSeq in (
            SELECT DISTINCT REF_MotifSeq
            FROM Instance
            WHERE REF_Motif=%s)
                ))
    """
    for REF_Motif in care.Motif.values():
        print REF_Motif
        care.cur.execute(SQL % (REF_Motif,REF_Motif))
    care.commit()

def usage():
    print "Stats.py -f <fasta_file> [-d] <dbfile> [-s]<seqname_list> [-m] <motif_list> [-p] pThresh [-o] <output>"

def has_file(filename):
    if not os.path.exists(filename):
        raise Exception,"There is no %s here" % filename

if __name__ == '__main__':
    (seqname_list,motif_list,pThresh)=(None,None,None)
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:s:m:p:o:", \
            ["help","database=","seqname_list=","motif_list=","pThresh","output="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-d", "--database"):
            dbname=arg
            has_file(dbname)
            output=dbname.split('.')[0]
        elif opt in ("-s", "--seqname_list"):
            seqname_list=arg
            has_file(seqname_list)
            output=seqname_list.split('.')[0]
            seqname_list=loadList(seqname_list)
        elif opt in ("-m", "--motif_list"):
            motif_list=arg
            has_file(motif_list)
            motif_list=loadList(motif_list)
        elif opt in ("-p", "--pThresh"):
            pThresh=float(arg)
        elif opt in ("-o", "--output"):
            output=arg
    if not 'dbname' in dir():
        print "tell me dbname please"
        usage()
        sys.exit(2)
    timestamp=time.time()
    #第一次执行的时候统计数据库中各类Motif对应的SeqName数量，以便看是否在给定列表中enrich
    Stats(dbname, seqname_list, motif_list, pThresh, output)
    print time.time()-timestamp
