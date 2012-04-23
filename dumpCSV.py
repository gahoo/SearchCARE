#! /usr/bin/env python
#coding=utf-8
import os
import sys
import time
import getopt
from zipfile import ZipFile, ZIP_DEFLATED
from CAREdb import CAREdb

def dumpCSV(dbname,SeqName_list):
    print "dumpCSV"
    CSV={}
    care=CAREdb(dbname)
    SeqName=invertDict(care.SeqName)
    Motif=invertDict(care.Motif)
    Organism=invertDict(care.Organism)

    if SeqName_list:
        REF_SeqName_list=[str(care.SeqName[sn]) for sn in SeqName_list if care.SeqName.has_key(sn)]
    else:
        REF_SeqName_list=None
    for rs in readSeqName(care,REF_SeqName_list):
        (REF_SeqName,REF_Motif,REF_Organism,Description,start,stop,strand,pValue)=rs
        if strand==1:
            strand='+'
        else:
            strand='-'
        Description=Description.replace(',',';')
        row=(Motif[REF_Motif],Organism[REF_Organism],Description,start,stop,strand,pValue)
        row=[str(c) for c in row]
        if CSV.has_key(SeqName[REF_SeqName]):
            CSV[SeqName[REF_SeqName]].append(",".join(row))
        else:
            CSV[SeqName[REF_SeqName]]=[",".join(row)]
    return CSV

def invertDict(d):
    return dict([(v,k) for k,v in d.items()])

def readSeqName(care,REF_SeqName_list):
    print "Load2RAM"
    care.cache('Instance')
    care.cache('Scanned',['REF_SeqName','REF_MotifSeq','start','stop','strand','pValue'],SeqName=REF_SeqName_list)
    if REF_SeqName_list:
        where="AND REF_SeqName in ('%s')" % "','".join(REF_SeqName_list)
    else:
        where=""
    SQL="""
    SELECT REF_SeqName,REF_Motif,REF_Organism,Description,start,stop,strand,pValue
    FROM cache.Instance,cache.Scanned
    WHERE cache.Instance.REF_MotifSeq=cache.Scanned.REF_MotifSeq
    %s
    """ % where
    #print SQL
    care.cur.execute(SQL)
    return care.cur.fetchall()

def writeCSV(outfile,CSV):
    print "writting CSV"
    zf=ZipFile(outfile, "w", ZIP_DEFLATED)
    for SeqName in CSV.keys():
        content="\n".join(CSV[SeqName])
        zf.writestr(SeqName+".csv",content)
    zf.close()

def usage():
    print "dumpCSV.py -d <dbfile> -i <SeqName_list> [-o] <outfile> "

def loadList(filename):
    return [line.strip() for line in open(filename).readlines()]

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:i:o:", ["help","database=", "filename=", "outfile="])
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
            SeqName_list=None
            outfile = "%s.zip" % dbname.split('.')[0]
        elif opt in ("-i", "--filename"):
            filename = arg
            if not os.path.exists(filename):
                raise Exception,"There is no %s here" % filename
            outfile = "%s.zip" % filename.split('.')[0]
            SeqName_list=loadList(filename)
        elif opt in ("-o", "--outfile"):
            outfile = arg
    if not 'dbname' in dir():
        usage()
        sys.exit(2)
    timestamp=time.time()
    CSV=dumpCSV(dbname,SeqName_list)
    writeCSV(outfile,CSV)
    print time.time()-timestamp
