#! /usr/bin/env python
#coding=utf-8
import os
import sys
import time
import getopt
from zipfile import ZipFile, ZIP_DEFLATED
from CAREdb import CAREdb
import cPickle

def dumpCoMotif(co_Occur,Merged_REF_SeqName_dist,output,zip=False):
    care=CAREdb(dbname)
    Seqnames=care.loadDB2dict('fa_file',0,1)
    if not os.path.exists(output):
        os.makedirs(output)
    if not os.path.exists("%s/gene_list" % output):
        os.makedirs("%s/gene_list" % output)
    #Co-occur
    for motifa,motifb in co_Occur.keys():
        print motifa,motifb
        filename="%s&%s.lst" % (motifa,motifb)
        outfile=open("%s/gene_list/%s" % (output,filename),'w')
        outfile.write("\n".join([Seqnames[k] for k in co_Occur[(motifa,motifb)].keys()]))
        outfile.close()
    #single motif
    motifs=[]
    [motifs.extend(k) for k in co_Occur.keys()]
    motifs=set(motifs)
    for motif in motifs:
        filename="%s.lst" % motif
        outfile=open("%s/gene_list/%s" % (output,filename),'w')
        outfile.write("\n".join([Seqnames[k] for k in Merged_REF_SeqName_dist[motif].keys()]))
        outfile.close()

def usage():
    print "dumpCoMotif.py -d <dbfile> -z <zpklfile> [-o] <outfile> -r 10#300"

def has_file(filename):
    if not os.path.exists(filename):
        raise Exception,"There is no %s here" % filename

if __name__ == '__main__':
    (output,co_dist_range)=(None,'10#300')
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:z:o:r:", ["help","database=", "zpklfile=", "output=","co_dist_range="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-d", "--database"):
            dbname=arg
            if not os.path.exists(dbname):
                raise Exception,"There is no %s here" % dbname
        elif opt in ("-z", "--zpklfile"):
            zpklfile=arg
            has_file(zpklfile)
            if not output:
                output=zpklfile.split('.')[0]
        elif opt in ("-o", "--output"):
            output=arg
        elif opt in ("-r", "--co_dist_range"):
            #arg=lower#higer
            co_dist_range=arg
    if not 'dbname' in dir() or not 'zpklfile' in dir():
        usage()
        sys.exit(2)

    timestamp=time.time()
    zf = ZipFile(zpklfile, 'r')
    Merged_REF_SeqName_dist=cPickle.loads(zf.open('Merged_REF_SeqName_dist.pkl').read())
    try:
        co_Occur=cPickle.loads(zf.open('co_Occur_%s.pkl' % co_dist_range).read())
    except KeyError:
        print "No co_Occur_%s.pkl in %s" % (co_dist_range,zpklfile)
        sys.exit(2)
    zf.close()
    dumpCoMotif(co_Occur,Merged_REF_SeqName_dist,output)
    print time.time()-timestamp
