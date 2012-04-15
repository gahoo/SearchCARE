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
from Stats import hypergeo_cdf,gauss_hypergeom,logchoose

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

def constrainDist(REF_SeqName,co_Occur,edge,motif_dist,co_dist_range,p=None,basename=None):
    '''
    co_Occur_calibrated[(motifa,motifb)][REF_SeqName]={(motifa_pos,motifb_pos):distance}
    '''
    (lower,higer)=co_dist_range
    edgefile="%s_coDist%d#%d.Edge" % (basename,lower,higer)
    EdgeOutput=open(edgefile, 'w')
    Seq_num=len(REF_SeqName)
    co_Occur_calibrated={}
    calibrated_edge=[]
    for row in edge:
        (motifa,motifb,co_num,a_num,b_num,hypergeocdf)=row
        co_Occur_calibrated[(motifa,motifb)]={}
        for REF_SeqName in co_Occur[(motifa,motifb)]:
            distA=motif_dist[motifa][REF_SeqName]
            distB=motif_dist[motifb][REF_SeqName]
            #print motifa,motifb,REF_SeqName
            distances=getCoDist(distA,distB,lower,higer)
            if distances:
            	print motifa,motifb,REF_SeqName,distances
                co_Occur_calibrated[(motifa,motifb)][REF_SeqName]=distances
        if not co_Occur_calibrated[(motifa,motifb)]:
            co_Occur_calibrated.pop((motifa,motifb))
            continue
        co_num=len(co_Occur_calibrated[(motifa,motifb)].keys())
        avg_dist=avgDist(co_Occur_calibrated[(motifa,motifb)])
        try:
            hypergeocdf=hypergeo_cdf(co_num,a_num,b_num,Seq_num)
        except AssertionError:
            print co_num,row
            continue
        if hypergeocdf<=p and co_num>=5:
            #row=(motifa,motifb,co_num,a_num,b_num,Seq_num,avg_dist,hypergeocdf)
            row=(motifa,motifb,co_num,a_num,b_num,hypergeocdf)
            calibrated_edge.append(row)
    EdgeOutput.write("\n".join(["\t".join([str(col) for col in row]) for row in calibrated_edge]))
    EdgeOutput.close()
    return co_Occur_calibrated

def usage():
    print "DistConstrain.py -z <zpkl> [-r] co_dist_range=lower#higer [-p] pThresh [-o] <output>"

def has_file(filename):
    if not os.path.exists(filename):
        raise Exception,"There is no %s here" % filename

if __name__ == '__main__':
    (pThresh,output,co_dist_range)=(0.05,None,(10,300))
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hz:r:p:o:", \
            ["help","zpklfile=","co_dist_range=","pThresh=","output="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-z", "--zpklfile"):
            zpklfile=arg
            has_file(zpklfile)
            if not output:
            	output=zpklfile.split('.')[0]
        elif opt in ("-r", "--co_dist_range"):
            #arg=lower#higer
            co_dist_range=map(int,arg.split("#"))
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
    constrainDist(REF_SeqName,co_Occur,edge,Merged_REF_SeqName_dist,co_dist_range,pThresh,output)
    zf.close()
    print time.time()-timestamp
