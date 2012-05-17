#! /usr/bin/env python
#coding=utf-8
import sys
import os
import time
import getopt
import cPickle
from zipfile import ZipFile, ZIP_DEFLATED
from Stats import loadDB2RAM,getREF_MotifSeq,getDist,getMerged,mergeREF_SeqName
from CAREdb import CAREdb

def getCoDist(distA,distB,lower,higer):
    '''
    获取符合距离范围的共现模体位置及距离
    '''
    comb=[(i,j) for i in distA for j in distB]
    distances=map(lambda x: x[0]-x[1], comb)
    CoDist={}
    for i,d in enumerate(distances):
        #print comb[i],d
        if abs(d)>=lower and abs(d)<=higer:
            CoDist[comb[i]]=d
    return CoDist

def getDBmotif_dist(dbname):
    care=CAREdb(dbname)
    print "Loading motif dist in DB..."
    (REF_SeqName,Motifs)=loadDB2RAM(care,None,None,False)
    id2MotifSeq=dict([(v,k) for k,v in care.MotifSeq.items()])
    motif_dist={}
    for i,motif in enumerate(Motifs):
        print motif,i
        REF_MotifSeq=getREF_MotifSeq(care,motif)
        motif_dist[motif]={}
        for MotifSeq_id in REF_MotifSeq:
            MotifSeq=id2MotifSeq[MotifSeq_id]
            motif_dist[motif][MotifSeq]=getDist(care,MotifSeq_id)
        print motif_dist[motif].keys()
    return getMerged(motif_dist,"REF_SeqName")

def getDBmotif_dist_zpkl(dbname):
    care=CAREdb(dbname)
    zpklfile="%s_motif_dist_db.zpkl" % dbname.split('.')[0]
    if os.path.exists(zpklfile):
        zf = ZipFile(zpklfile, 'r', ZIP_DEFLATED)
        motifs_zipped=[m.split('.')[0] for m in zf.namelist()]
        zf.close()
        print motifs_zipped
    else:
        motifs_zipped=[]
    zf = ZipFile(zpklfile, 'a', ZIP_DEFLATED)
    print "Loading motif dist in DB..."
    (REF_SeqName,Motifs)=loadDB2RAM(care,None,None,False)
    id2MotifSeq=dict([(v,k) for k,v in care.MotifSeq.items()])
    for i,motif in enumerate(Motifs):
        print motif,i
        if motif in motifs_zipped:
            continue
        REF_MotifSeq=getREF_MotifSeq(care,motif)
        motif_dist={}
        for MotifSeq_id in REF_MotifSeq:
            MotifSeq=id2MotifSeq[MotifSeq_id]
            motif_dist[MotifSeq]=getDist(care,MotifSeq_id)
        motif_dist=mergeREF_SeqName(motif_dist)
        print motif_dist.keys()
        zf.writestr('%s.pkl' % motif, cPickle.dumps(motif_dist))
    zf.close()
    motif_dist=loadZpkl_motif_dist(zpklfile,Motifs)
    return motif_dist #getMerged(motif_dist,"REF_SeqName")

def loadZpkl_motif_dist(zpklfile,Motifs):
    zf = ZipFile(zpklfile, 'r')
    motif_dist={}
    for i,motif in enumerate(Motifs):
        print motif,i
        motif_dist[motif] = cPickle.loads(zf.open('%s.pkl' % motif).read())
    return motif_dist

def getCoOccurMotifDB(motif_dist_db,co_dist_range):
    Motifs=motif_dist_db.keys()
    co_Occur_db={}
    for motifa in Motifs:
        print len(Motifs)
        for motifb in Motifs[1:]:
            #print motifa,motifb
            distances=getCoMotifDist(motifa,motifb,motif_dist_db,co_dist_range)
            if distances:
                #co_Occur_db[(motifa,motifb)]=distances
                #为了减小内存占用，仅保存数量信息
                co_Occur_db[(motifa,motifb)]=len(distances.keys())
        Motifs=Motifs[1:]
    return co_Occur_db

def getCoMotifDist(motifa,motifb,motif_dist_db,co_dist_range):
    (lower,higer)=co_dist_range
    REF_SeqNames=set(motif_dist_db[motifa].keys())&set(motif_dist_db[motifb].keys())
    motif_pair_dist={}
    for REF_SeqName in REF_SeqNames:
        distA_db=motif_dist_db[motifa][REF_SeqName]
        distB_db=motif_dist_db[motifb][REF_SeqName]
        distances_db=getCoDist(distA_db,distB_db,lower,higer)
        if distances_db:
            motif_pair_dist[REF_SeqName]=distances_db
    return motif_pair_dist

def getCoOccurMotifDBZpkl(zpklfile,co_dist_range):
    zf = ZipFile(zpklfile, 'r', ZIP_DEFLATED)
    Motifs=[m.split('.')[0] for m in zf.namelist()]
    co_Occur_db={}
    for motifa in Motifs:
        print len(Motifs)
        motif_dist_db_motifa=cPickle.loads(zf.open('%s.pkl' % motifa).read())
        for motifb in Motifs[1:]:
            #print motifa,motifb
            motif_dist_db_motifb=cPickle.loads(zf.open('%s.pkl' % motifb).read())
            distances=getCoMotifDist(motif_dist_db_motifa,motif_dist_db_motifb,co_dist_range)
            if distances:
                #co_Occur_db[(motifa,motifb)]=distances
                #为了减小内存占用，仅保存数量信息
                co_Occur_db[(motifa,motifb)]=len(distances.keys())
        Motifs=Motifs[1:]
    zf.close()
    return co_Occur_db

def getCoMotifDistZpkl(motif_dist_db_motifa,motif_dist_db_motifb,co_dist_range):
    (lower,higer)=co_dist_range
    REF_SeqNames=set(motif_dist_db_motifa.keys())&set(motif_dist_db_motifb.keys())
    motif_pair_dist={}
    for REF_SeqName in REF_SeqNames:
        distA_db=motif_dist_db_motifa[REF_SeqName]
        distB_db=motif_dist_db_motifb[REF_SeqName]
        distances_db=getCoDist(distA_db,distB_db,lower,higer)
        if distances_db:
            motif_pair_dist[REF_SeqName]=distances_db
    return motif_pair_dist

def slimCoOccurMotifDB(co_Occur_db):
    co_Occur_slim={}
    for motif_pair in co_Occur_db.keys():
        co_Occur_slim[motif_pair]=len(co_Occur_db[motif_pair].keys())
    return co_Occur_slim

def savePKL(pklname,data):
    pkl=open(pklname,'w')
    cPickle.dump(data,pkl)
    pkl.close()

def usage():
    print "dumpDBmotifDist.py -d <dbname> -r 10#250"

def has_file(filename):
    if not os.path.exists(filename):
        raise Exception,"There is no %s here" % filename

if __name__ == '__main__':
    co_dist_range=(10,250)
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:r:", \
            ["help","database="])
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
        elif opt in ("-r", "--co_dist_range"):
            #arg=lower#higer
            co_dist_range=tuple(map(int,arg.split("#")))
    if not 'dbname' in dir():
        print "tell me dbname please"
        usage()
        sys.exit(2)

    timestamp=time.time()

    motif_dist_pkl='%s_motif_dist_db.pkl' % output
    if os.path.exists(motif_dist_pkl):
        motif_dist_db=cPickle.loads(open(motif_dist_pkl).read())
    else:
        motif_dist_db=getDBmotif_dist_zpkl(dbname)
        savePKL(motif_dist_pkl,motif_dist_db)

    co_Occur_pkl='%s_co_Occur_db%s#%s.pkl' % (output,co_dist_range[0],co_dist_range[1])
    co_Occur_db=getCoOccurMotifDB(motif_dist_db,co_dist_range)
    savePKL(co_Occur_pkl,co_Occur_db)
    '''
    co_Occur_slim_pkl='%s_co_Occur_slim%s#%s.pkl' % (output,co_dist_range[0],co_dist_range[1])
    slimCoOccurMotifDB=slimCoOccurMotifDB(co_Occur_db)
    savePKL(co_Occur_slim_pkl,slimCoOccurMotifDB)
    '''
    print time.time()-timestamp