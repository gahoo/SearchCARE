#! /usr/bin/env python
#coding=utf-8
import sys
import os
import time
from CAREdb import CAREdb
from CAREdb import loadMotifs2list
from subprocess import Popen,PIPE
from zipfile import ZipFile, ZIP_DEFLATED
import matplotlib
matplotlib.use("Agg")
import pylab
import cPickle

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

def loadScanned(scanned_file):
    '''
    加载扫描完的文件
    '''
    zf = ZipFile(scanned_file, 'r')
    scanned_file=zf.namelist()[0]
    scanned = [rs.split('\t') for rs in zf.read(scanned_file).split('\n')[1:-1]]
    zf.close()
    return scanned

def SearchCARE(input='test.fasta'):
    '''
    在序列中搜索顺式调控元件motif
    '''
    basename=input.split('.')[0]
    CARE=CAREdb('test.db')
    Motif_lst=loadMotifs2list('test.db')
    scanned_path="%s/Scan" % basename
    if not os.path.exists(scanned_path):
        os.makedirs(scanned_path)
    scanned=checkScanned("%s/Scan/" % basename)
    for entry in Motif_lst:
        (motif_name, Description, Sequence)=entry
        Description=Description.replace(",",";")
        print motif_name,Sequence
        if Sequence in scanned:
            print "Skipping"
            continue
        ScanMotif(Sequence,input)
    parseScan(input=input)

def checkScanned(chk_path):
    '''
    检查已经扫描完毕的文件列表
    返回已经扫描完毕的序列
    '''
    return [rs.split('.')[0] for rs in os.listdir(chk_path)]

def parseScan(input='test.fasta', zip=True):
    '''
    dump[motif_name][Sequence][LOC][strand]=[start,end]
    '''
    basename=input.split('.')[0]
    Motif_lst=loadMotifs2list('test.db')
    CSV={}
    Motif_dic={}
    for entry in Motif_lst:
        (motif_name,Sequence)=(entry[0],entry[2])
        result=loadScanned("%s/Scan/%s.zip" % (basename,Sequence))
        seq_dict=ScanParser(entry,result,CSV)
        #TODO
        nestDictUpdate(Motif_dic,[motif_name,Sequence],seq_dict)
    if zip:
        file=ZipFile("%s_csv.zip" % basename, "w", ZIP_DEFLATED)
        for LOC in CSV.keys():
            file.writestr(LOC+".csv","\n".join(CSV[LOC]))
        file.close()
    else:
        for LOC in CSV.keys():
            writeCSV(LOC,"\n".join(CSV[LOC]))
    zf = ZipFile("%s.zpkl" % basename, 'w', ZIP_DEFLATED)
    zf.writestr('Motif_dic.pkl', cPickle.dumps(Motif_dic))
    zf.close()

def ScanParser(entry,scan_res,csv):
    '''
    解析FIMO结果，得到LOC列表，构建CSV，motif分布位置信息
    seq_dict[LOC][strand]=[(start,end),]
    '''
    (motif_name, Description, Sequence)=entry
    Description=Description.replace(",",";")
    print motif_name,Sequence
    seq_dict={}
    for rs in scan_res:
        (LOC, start, stop, strand, Sequence) = \
        (rs[1],int(rs[2]),int(rs[3]),rs[4],rs[7])
        length=str(abs(start-stop))
        row=",".join([motif_name,str(start),strand,length,Sequence,Description])
        nestDictList(csv,[LOC],row)
        nestDictList(seq_dict,[LOC,strand],(start,stop))
    return seq_dict

def nestDictList(dic,keys,value):
    '''
    递归设置字典，即直接赋值dic[k1][k2][k3]=v，当k不存在时，自动创建
    nestDict(seq_dict,[LOC,strand],(start,stop))
    相当于seq_dict[LOC][strand].append(start,stop)
    '''
    if len(keys)==1:
        if not dic.has_key(keys[0]):
            dic[keys[0]]=[value]
        else:
            dic[keys[0]].append(value)
        return
    if not dic.has_key(keys[0]):
        dic[keys[0]]={}
    nestDictList(dic[keys[0]],keys[1:],value)

def nestDictUpdate(dic,keys,value):
    '''
    递归更新字典，即直接赋值dic[k1][k2][k3]=v，当k不存在时，自动创建
    nestDictUpdate(Motif_dic,[motif_name,Sequence],seq_dict)
    相当于Motif_dic[motif_name][Sequence].update(seq_dict)
    '''
    if len(keys)==1:
        if not dic.has_key(keys[0]):
            dic[keys[0]]=value
        else:
            dic[keys[0]].update(value)
        return
    if not dic.has_key(keys[0]):
        dic[keys[0]]={}
    nestDictUpdate(dic[keys[0]],keys[1:],value)

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
    Motif_LOC=SearchCARE(input=sys.argv[1])
    print time.time()-timestamp
