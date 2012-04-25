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
from Stats import hypergeo_cdf,gauss_hypergeom,logchoose,loadDB2RAM,getREF_MotifSeq,getDist,getMerged,hypergeo_cdf_PAN
from CAREdb import CAREdb
from dumpDBmotifDist import getDBmotif_dist

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

def avgDist(co_dict):
    distances=[v.values()[0] for v in co_dict.values()]
    distances=map(abs,distances)
    avg_dist=sum(distances)/len(distances)
    return avg_dist

'''
def getDBmotif_dist(care,edge):
    print "Loading motif dist in DB..."
    Motifs=[]
    for row in edge:
        Motifs.extend(row[:2])
    Motifs=set(Motifs)
    (REF_SeqName,Motifs)=loadDB2RAM(care,None,Motifs,False)
    id2MotifSeq=dict([(v,k) for k,v in care.MotifSeq.items()])
    motif_dist={}
    for motif in Motifs:
        print motif
        REF_MotifSeq=getREF_MotifSeq(care,motif)
        motif_dist[motif]={}
        for MotifSeq_id in REF_MotifSeq:
            MotifSeq=id2MotifSeq[MotifSeq_id]
            motif_dist[motif][MotifSeq]=getDist(care,MotifSeq_id)
    return getMerged(motif_dist,"REF_SeqName")
'''


def constrainDistPAN(dbname,motif_dist_dumps,REF_SeqName,co_Occur,edge,motif_dist,co_dist_range,p=None,basename=None):
    '''
    co_Occur_calibrated[(motifa,motifb)][REF_SeqName]={(motifa_pos,motifb_pos):distance}
    '''
    if motif_dist_dumps:
        print "load motif_dist_dumps from PKL dumps"
        motif_dist_db=cPickle.loads(open(motif_dist_dumps).read())
    else:
        print "load motif_dist_dumps from db"
        motif_dist_db=getDBmotif_dist(dbname)
    (lower,higer)=co_dist_range
    list_size=len(REF_SeqName)
    care=CAREdb(dbname)
    db_size=len(care.SeqName)
    co_Occur_calibrated={}
    co_Occur_calibrated_db={}
    calibrated_edge=[]
    for row in edge:
        (motifa,motifb,co_num,a_num,b_num,hypergeocdf)=row
        co_Occur_calibrated[(motifa,motifb)]={}
        co_Occur_calibrated_db[(motifa,motifb)]={}
        #print motifa,motifb
        for REF_SeqName in co_Occur[(motifa,motifb)]:
            distA=motif_dist[motifa][REF_SeqName]
            distB=motif_dist[motifb][REF_SeqName]
            #print motifa,motifb,REF_SeqName
            distances=getCoDist(distA,distB,lower,higer)
            if distances:
                #print motifa,motifb,REF_SeqName,distances
                co_Occur_calibrated[(motifa,motifb)][REF_SeqName]=distances
        if not co_Occur_calibrated[(motifa,motifb)]:
            co_Occur_calibrated.pop((motifa,motifb))
            co_Occur_calibrated_db.pop((motifa,motifb))
            continue

        REF_SeqNames=set(motif_dist_db[motifa].keys())&set(motif_dist_db[motifb].keys())
        for REF_SeqName in REF_SeqNames:
            distA_db=motif_dist_db[motifa][REF_SeqName]
            distB_db=motif_dist_db[motifb][REF_SeqName]
            distances_db=getCoDist(distA_db,distB_db,lower,higer)
            if distances_db:
                co_Occur_calibrated_db[(motifa,motifb)][REF_SeqName]=distances_db
        co_num=len(co_Occur_calibrated[(motifa,motifb)].keys())
        co_num_db=len(co_Occur_calibrated_db[(motifa,motifb)].keys())
        print co_num_db,len(REF_SeqNames)
        avg_dist=avgDist(co_Occur_calibrated[(motifa,motifb)])
        try:
            hypergeocdf=hypergeo_cdf_PAN(co_num,list_size,co_num_db,db_size)
        except AssertionError:
            print co_num,list_size,co_num_db,db_size
            co_Occur_calibrated.pop((motifa,motifb))
            co_Occur_calibrated_db.pop((motifa,motifb))
            continue
        if hypergeocdf<=p and co_num>=5:
            #row=(motifa,motifb,co_num,a_num,b_num,list_size,avg_dist,hypergeocdf)
            row=(motifa,motifb,co_num,list_size,co_num_db,hypergeocdf)
            calibrated_edge.append(row)
        else:
            print (motifa,motifb,co_num,list_size,co_num_db,hypergeocdf)
            co_Occur_calibrated.pop((motifa,motifb))

    if calibrated_edge:
        edgefile="%s_coDist%d#%d.Edge" % (basename,lower,higer)
        EdgeOutput=open(edgefile, 'w')
        EdgeOutput.write("\n".join(["\t".join([str(col) for col in row]) for row in calibrated_edge]))
        EdgeOutput.close()
    else:
        print "no significant edges"
        return
    return co_Occur_calibrated

def constrainDist(REF_SeqName,co_Occur,edge,motif_dist,co_dist_range,p=None,basename=None):
    '''
    co_Occur_calibrated[(motifa,motifb)][REF_SeqName]={(motifa_pos,motifb_pos):distance}
    '''
    (lower,higer)=co_dist_range
    Seq_num=len(REF_SeqName)
    co_Occur_calibrated={}
    calibrated_edge=[]
    for row in edge:
        (motifa,motifb,co_num,a_num,b_num,hypergeocdf)=row
        co_Occur_calibrated[(motifa,motifb)]={}
        #print motifa,motifb
        for REF_SeqName in co_Occur[(motifa,motifb)]:
            distA=motif_dist[motifa][REF_SeqName]
            distB=motif_dist[motifb][REF_SeqName]
            #print motifa,motifb,REF_SeqName
            distances=getCoDist(distA,distB,lower,higer)
            if distances:
                #print motifa,motifb,REF_SeqName,distances
                co_Occur_calibrated[(motifa,motifb)][REF_SeqName]=distances
        if not co_Occur_calibrated[(motifa,motifb)]:
            co_Occur_calibrated.pop((motifa,motifb))
            continue
        co_num=len(co_Occur_calibrated[(motifa,motifb)].keys())
        avg_dist=avgDist(co_Occur_calibrated[(motifa,motifb)])
        try:
            hypergeocdf=hypergeo_cdf(co_num,a_num,b_num,Seq_num)
        except AssertionError:
            print co_num,row,Seq_num
            co_Occur_calibrated.pop((motifa,motifb))
            continue
        if hypergeocdf<=p and co_num>=5:
            #row=(motifa,motifb,co_num,a_num,b_num,Seq_num,avg_dist,hypergeocdf)
            row=(motifa,motifb,co_num,a_num,b_num,hypergeocdf)
            calibrated_edge.append(row)
        else:
            print (motifa,motifb,co_num,a_num,b_num,hypergeocdf)
            co_Occur_calibrated.pop((motifa,motifb))

    if calibrated_edge:
        edgefile="%s_coDist%d#%d.Edge" % (basename,lower,higer)
        EdgeOutput=open(edgefile, 'w')
        EdgeOutput.write("\n".join(["\t".join([str(col) for col in row]) for row in calibrated_edge]))
        EdgeOutput.close()
    else:
        print "no significant edges"
        return
    return co_Occur_calibrated

def usage():
    print "DistConstrain.py -z <zpkl> [-r] co_dist_range=lower#higer [-p] pThresh [-o] <output>"

def has_file(filename):
    if not os.path.exists(filename):
        raise Exception,"There is no %s here" % filename

if __name__ == '__main__':
    (pThresh,output,co_dist_range,motif_dist_dumps)=(0.05,None,(10,250),None)
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:m:z:r:p:o:", \
            ["help","database=","motif_dist_dumps=","zpklfile=","co_dist_range=","pThresh=","output="])
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
        elif opt in ("-m", "--motif_dist_dumps"):
            motif_dist_dumps=arg
            has_file(motif_dist_dumps)
        elif opt in ("-z", "--zpklfile"):
            zpklfile=arg
            has_file(zpklfile)
            if not output:
                output=zpklfile.split('.')[0]
        elif opt in ("-r", "--co_dist_range"):
            #arg=lower#higer
            co_dist_range=tuple(map(int,arg.split("#")))
        elif opt in ("-p", "--pThresh"):
            pThresh=float(arg)
        elif opt in ("-o", "--output"):
            output=arg
    if not 'zpklfile' in dir():
        print "tell me zpklfile please"
        usage()
        sys.exit(2)
    
    timestamp=time.time()
    zf = ZipFile(zpklfile, 'r')
    REF_SeqName = cPickle.loads(zf.open('REF_SeqName.pkl').read())
    co_Occur = cPickle.loads(zf.open('co_Occur.pkl').read())
    edge=cPickle.loads(zf.open('edge.pkl').read())
    Merged_REF_SeqName_dist=cPickle.loads(zf.open('Merged_REF_SeqName_dist.pkl').read())
    co_Occur_calibrated=constrainDistPAN(dbname,motif_dist_dumps,REF_SeqName,co_Occur,edge,Merged_REF_SeqName_dist,co_dist_range,pThresh,output)
    #co_Occur_calibrated=constrainDist(REF_SeqName,co_Occur,edge,Merged_REF_SeqName_dist,co_dist_range,pThresh,output)
    zf.close()
    zf = ZipFile(zpklfile, 'a')
    if co_Occur_calibrated:
        zf.writestr('co_Occur_%s#%s.pkl' % co_dist_range , cPickle.dumps(co_Occur_calibrated))
    zf.close()
    print time.time()-timestamp
