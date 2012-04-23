#! /usr/bin/env python
#coding=utf-8
import sys
import os
import time
import getopt
from CAREdb import CAREdb
from subprocess import Popen,PIPE
from zipfile import ZipFile, ZIP_DEFLATED

def ScanMotif(Motif, fasta_file, bgfile):
    '''
    Motif:iupac格式的字符串
    fasta_file:待寻找的fasta文件
    '''
    iupac2meme=Popen(args=['iupac2meme',Motif],stdout=PIPE)
    (motif,err)=iupac2meme.communicate()
    #iupac2meme.wait()
    fimo_arg=['fimo']
    if bgfile:
        fimo_arg.extend(['--bgfile', bgfile])
    fimo_arg.extend(['--no-qvalue','-text','--output-pthresh','8e-5','--verbosity','1','-', fasta_file])
    fimo=Popen(args=fimo_arg,stdin=PIPE,stdout=PIPE)
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
    SELECT DISTINCT Motif,Sequence
    FROM Motif,Instance,MotifSeq
    WHERE Motif.id=Instance.REF_Motif
        AND MotifSeq.id=Instance.REF_MotifSeq
        AND MotifSeq.id NOT IN ( '%s' )
    """ % Skip
    #print SQL
    caredb.cur.execute(SQL)
    return caredb.cur.fetchall()

def SearchCARE(filename,dbname,bgfile):
    '''
    在序列中搜索顺式调控元件motif
    '''
    basename = filename.split('.')[0]
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
        CARE.addScanned(ScanMotif(Sequence,filename,bgfile))
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

def usage():
    print "buildDB.py -i <fasta_file> [-d] <dbfile> [-b]<bgfile>"

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:d:b:", ["help", "filename=", "database=", "bgfile="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--filename"):
            filename = arg
            if not os.path.exists(filename):
                raise Exception,"There is no %s here" % filename
            dbname = "%s.db" % filename.split('.')[0]
            bgfile=None
        elif opt in ("-d", "--database"):
            dbname=arg
        elif opt in ("-b", "--bgfile"):
            bgfile = arg
    if not 'dbname' in dir():
        usage()
        sys.exit(2)
    timestamp=time.time()
    SearchCARE(filename,dbname,bgfile)
    print time.time()-timestamp