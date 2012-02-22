#! /usr/bin/env python
#coding=utf-8
import sys
import os
import time
from CAREdb import CAREdb
from subprocess import Popen,PIPE
from zipfile import ZipFile, ZIP_DEFLATED

def ScanMotif(Motif, fasta_file):
    '''
    Motif:iupac格式的字符串
    fasta_file:待寻找的fasta文件
    '''
    iupac2meme=Popen(args=['iupac2meme',Motif],stdout=PIPE)
    (motif,err)=iupac2meme.communicate()
    #iupac2meme.wait()
    fimo=Popen(args=['fimo','--no-qvalue','-text','--output-pthresh','8e-5','--verbosity','1','-', fasta_file],stdin=PIPE,stdout=PIPE)
    (result,err)=fimo.communicate(motif)
    #fimo.wait()
    basename=fasta_file.split('.')[0]
    dumpzfile="%s/Scan/%s.zip" % (basename,Motif)
    zf=ZipFile(dumpzfile, "w", ZIP_DEFLATED)
    zf.writestr(Motif+".scan",result)
    zf.close()
    return [rs.split('\t') for rs in result.split('\n')[1:-1]]

def loadMotifs(caredb):
    '''
    加载还没有扫描过的Motif序列
    '''
    Skip="','".join(checkScannedDB(caredb))
    SQL="""
    SELECT Motif,Sequence
    FROM Motif,Instance,MotifSeq
    WHERE Motif.id=Instance.REF_Motif
        AND MotifSeq.id=Instance.REF_MotifSeq
        AND MotifSeq.id NOT IN ( '%s' )
    """ % Skip
    #print SQL
    caredb.cur.execute(SQL)
    return caredb.cur.fetchall()

def SearchCARE(filename):
    '''
    在序列中搜索顺式调控元件motif
    '''
    basename=filename.split('.')[0]
    dbname='%s.db' % basename
    if not os.path.exists(dbname):
        CARE=CAREdb(dbname)
        CARE.importMotifs('CARE.txt')
    else:
        CARE=CAREdb(dbname)
    Motif_lst=loadMotifs(CARE)
    scanned=checkScannedLog("%s.log" % basename)
    scanlog=open("%s.log" % basename,'a')
    scanned_path="%s/Scan" % basename
    if not os.path.exists(scanned_path):
        os.makedirs(scanned_path)
    for entry in Motif_lst:
        (motif_name, Sequence)=entry
        print motif_name,Sequence
        if Sequence in scanned:
            print "skip"
            continue
        CARE.addScanned(ScanMotif(Sequence,filename))
        scanned.append(Sequence)
        scanlog.write(Sequence+'\n')

def checkScannedDB(caredb):
    '''
    检查已经扫描完毕的文件列表
    返回已经扫描完毕的序列
    '''
    SQL="""
    SELECT distinct REF_MotifSeq
    FROM scanned
    """
    caredb.cur.execute(SQL)
    return [str(rs[0]) for rs in caredb.cur.fetchall()]
    #return [rs.split('.')[0] for rs in os.listdir(chk_path)]

def checkScannedLog(logfile):
    '''
    检查并加载已经扫描完毕的序列
    '''
    if os.path.exists(logfile):
        return open(logfile,'r').read().split('\n')
    else:
        return []

if __name__ == "__main__":
    try:
        sys.argv[1]
    except:
        print "Tell me the filename please."
        sys.exit()
    if not os.path.exists(sys.argv[1]):
        raise Exception,"There is no file name %s" % sys.argv[1]
    basename=sys.argv[1].split('.')[0]
    timestamp=time.time()
    SearchCARE(sys.argv[1])
    print time.time()-timestamp
