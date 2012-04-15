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

def Stats(dbname, SeqName=None, Motifs=None, co_dist_range=None, p=None, output=None):
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
    Merged_MotifSeq_dist=getMerged(motif_dist,'MotifSeq')
    Merged_REF_SeqName_dist=getMerged(motif_dist,"REF_SeqName")
    SeqName_counts=countSeqName(MergedSeqName)
    zf = ZipFile("%s.zpkl" % basename, 'w', ZIP_DEFLATED)
    if SeqName:
        Enriched=enrichMent(care,SeqName_counts,list_size=len(REF_SeqName),output=output)
        zf.writestr('Enriched.pkl', cPickle.dumps(Enriched))
    zf.writestr('REF_SeqName.pkl', cPickle.dumps(REF_SeqName))
    zf.writestr('motif_dist.pkl', cPickle.dumps(motif_dist))
    zf.writestr('Merged_MotifSeq_dist.pkl', cPickle.dumps(Merged_MotifSeq_dist))
    zf.writestr('Merged_dist.pkl', cPickle.dumps(Merged_dist))
    zf.writestr('Motif_SeqName.pkl', cPickle.dumps(Motif_SeqName))
    zf.writestr('MergedSeqName.pkl', cPickle.dumps(MergedSeqName))
    zf.writestr('Merged_REF_SeqName_dist.pkl', cPickle.dumps(Merged_REF_SeqName_dist))
    zf.writestr('SeqName_counts.pkl', cPickle.dumps(SeqName_counts))
    (co_Occur,edge)=coOccur(care,REF_SeqName,SeqName_counts,MergedSeqName,p,output)
    zf.writestr('co_Occur.pkl', cPickle.dumps(co_Occur))
    zf.writestr('edge.pkl', cPickle.dumps(edge))
    if co_dist_range:
        co_Occur_calibrated=calibrateDist(care,REF_SeqName,co_Occur,edge,motif_dist,co_dist_range,p,output)
        zf.writestr('co_Occur_calibrated.pkl', cPickle.dumps(co_Occur_calibrated))
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
    dist[REF_SeqName][strand]=start/stop
    '''
    if strand=='both':
        WHERE=""
    elif strand=='+':
        WHERE="AND strand=1"
    elif strand=='-':
        WHERE="AND strand=0"
    SQL="""
    SELECT REF_SeqName,start,stop,strand
    FROM cache.Scanned
    WHERE REF_MotifSeq=%s
        %s
    """ % (MotifSeq_id,WHERE)
    care.cur.execute(SQL)
    rs=care.cur.fetchall()
    dist=dict([(r[0],{'+':[],'-':[]}) for r in rs])
    for REF_SeqName,start,stop,strand in rs:
        if strand==1:
            dist[REF_SeqName]['+'].append(start)
        else:
            dist[REF_SeqName]['-'].append(stop)
    return dist

def getCoDist(distA,distB,lower,higer):
    '''
    获取符合距离范围的共现模体位置及距离
    '''
    comb=[(i,j) for i in distA for j in distB]
    distances=map(lambda x: x[0]-x[1], comb)
    CoDist={}
    for i,d in enumerate(distances):
        print comb[i],d
        if abs(d)>=lower and abs(d)<=higer:
            CoDist[comb[i]]=d
    return CoDist

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
        #print x,m,list_size,db_size,ER,stats.hypergeom.cdf(x,db_size,m,list_size)
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
        if p<0.05 and x>=5:
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
            co_Occur[(motifa,motifb)]=motifA&motifB
            co_num=len(motifA&motifB)
            a_num=SeqName_counts[motifa]
            b_num=SeqName_counts[motifb]
            hypergeocdf=hypergeo_cdf(co_num,a_num,b_num,Seq_num)
            #共现小于5已经没有任何意义，应该用Fisher精确检验
            if hypergeocdf<=p and co_num>=5:
                row=(motifa,motifb,co_num,a_num,b_num,hypergeocdf)
                #print row
                edge.append(row)
    EdgeOutput.write("\n".join(["\t".join([str(col) for col in row]) for row in edge]))
    EdgeOutput.close()
    #EdgeOutput=open("%s_cali.Edge" % basename, 'w')
    return (co_Occur,edge)

def calibrateDist(care,REF_SeqName,co_Occur,edge,motif_dist,co_dist_range,p=None,output=None):
    '''
    co_Occur_calibrated[(motifa,motifb)][REF_SeqName]={(motifa_pos,motifb_pos):distance}
    '''
    (lower,higer)=co_dist_range
    if not p:
        p=0.05
    if not output:
        basename=care.dbfile.split('.')[0]
    else:
        basename=output
    edgefile="%s_codist%d#%d.Edge" % (basename,lower,higer)
    EdgeOutput=open(edgefile, 'w')
    if REF_SeqName:
        Seq_num=len(REF_SeqName)
    else:
        Seq_num=len(care.SeqName)
    dist=getMerged(motif_dist,"REF_SeqName")
    co_Occur_calibrated={}
    calibrated_edge=[]
    for row in edge:
        (motifa,motifb,co_num,a_num,b_num,hypergeocdf)=row
        co_Occur_calibrated[(motifa,motifb)]={}
        for REF_SeqName in co_Occur[(motifa,motifb)]:
            distA=dist[motifa][REF_SeqName]
            distB=dist[motifb][REF_SeqName]
            print motifa,motifb,REF_SeqName
            distances=getCoDist(distA,distB,lower,higer)
            if distances:
                co_Occur_calibrated[(motifa,motifb)][REF_SeqName]=distances
        if not co_Occur_calibrated[(motifa,motifb)]:
            co_Occur_calibrated.pop((motifa,motifb))
            continue
        co_num=len(co_Occur_calibrated[(motifa,motifb)].keys())
        avg_dist=avgDist(co_Occur_calibrated[(motifa,motifb)])
        try:
            hypergeocdf=hypergeo_cdf(co_num,a_num,b_num,Seq_num)
        except:
            print co_num,row
            continue
        if hypergeocdf<=p and co_num>=5:
            #row=(motifa,motifb,co_num,a_num,b_num,Seq_num,avg_dist,hypergeocdf)
            row=(motifa,motifb,co_num,a_num,b_num,hypergeocdf)
            calibrated_edge.append(row)
    EdgeOutput.write("\n".join(["\t".join([str(col) for col in row]) for row in calibrated_edge]))
    return co_Occur_calibrated

def getMerged(motif_dict,type):
    '''
    合并Motif_SeqName结果，即把：
    Motif_SeqName[motif][MotifSeq]={'+':SeqName,'-':SeqName}
    变成：
    Motif_SeqName[motif]=SeqName
    或
    motif_dist[motif][MotifSeq][REF_SeqName]={'+':start,'-':stop}
    变成：
    motif_dist[motif]=dist
    或
    motif_dist[motif][MotifSeq][REF_SeqName]={'+':start,'-':stop}
    变成：
    motif_dist[motif][MotifSeq]={'+':start,'-':stop}
    或
    motif_dist[motif][MotifSeq][REF_SeqName]={'+':start,'-':stop}
    变成：
    motif_dist[motif][REF_SeqName]=dist
    '''
    Merged={}
    if type=='dist':
        for motif in motif_dict.keys():
            Merged[motif]=[]
            for MotifSeq in motif_dict[motif].keys():
                Merged[motif].extend(merge(motif_dict[motif][MotifSeq]))
    elif type=='SeqName':
        for motif in motif_dict.keys():
            Merged[motif]=mergeSeqName(motif_dict[motif])
    elif type=='REF_SeqName':
        for motif in motif_dict.keys():
            Merged[motif]=mergeREF_SeqName(motif_dict[motif])
    elif type=='MotifSeq':
        for motif in motif_dict.keys():
            Merged[motif]={}
            for MotifSeq in motif_dict[motif].keys():
                Merged[motif][MotifSeq]=mergeMotifSeq(motif_dict[motif][MotifSeq])
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

def mergeREF_SeqName(motif_dict):
    '''
    motif[MotifSeq][REF_SeqName]={'+':start,'-':stop}
    return motif[REF_SeqName]=dist
    '''
    REF_SeqName_merged={}
    REF_SeqName_keys=[]
    for MotifSeq in motif_dict.keys():
        REF_SeqName_keys.extend(motif_dict[MotifSeq].keys())
    for REF_SeqName in set(REF_SeqName_keys):
        REF_SeqName_merged[REF_SeqName]=[]
    for MotifSeq in motif_dict.keys():
        for REF_SeqName in motif_dict[MotifSeq].keys():
            d=motif_dict[MotifSeq][REF_SeqName]
            m=[]
            m.extend(d['+']+d['-'])
            REF_SeqName_merged[REF_SeqName].extend(m)
    return REF_SeqName_merged

def mergeMotifSeq(motif_dict):
    '''
    motif_dict[REF_SeqName]={'+':start,'-':stop}
    motif_dict={'+':start,'-':stop}
    '''
    MotifSeq_merge={'+':[],'-':[]}
    for REF_SeqName in motif_dict.keys():
        if motif_dict[REF_SeqName]:
            MotifSeq_merge['+'].extend(motif_dict[REF_SeqName]['+'])
            MotifSeq_merge['-'].extend(motif_dict[REF_SeqName]['-'])
    return MotifSeq_merge

def countSeqName(MergedSeqName):
    '''
    统计SeqName数目
    '''
    SeqName_counts={}
    for motif in MergedSeqName.keys():
        SeqName_counts[motif]=len(MergedSeqName[motif])
    return SeqName_counts

def avgDist(co_dict):
    distances=[v.values()[0] for v in co_dict.values()]
    distances=map(abs,distances)
    avg_dist=sum(distances)/len(distances)
    return avg_dist

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
    print "Stats.py -f <fasta_file> [-d] <dbfile> [-s]<seqname_list> [-m] <motif_list> [-r] co_dist_range=lower#higer [-p] pThresh [-o] <output>"

def has_file(filename):
    if not os.path.exists(filename):
        raise Exception,"There is no %s here" % filename

if __name__ == '__main__':
    (seqname_list,motif_list,pThresh,co_dist_range)=(None,None,None,(10,300))
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:s:m:r:p:o:", \
            ["help","database=","seqname_list=","motif_list=","co_dist_range=","pThresh=","output="])
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
        elif opt in ("-r", "--co_dist_range"):
            #arg=lower#higer
            co_dist_range=map(int,arg.split("#"))
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
    Stats(dbname, seqname_list, motif_list, co_dist_range, pThresh, output)
    print time.time()-timestamp
