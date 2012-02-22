#! /usr/bin/env python
#coding=utf-8
import os
import sqlite3
import sys
from sqlitedb import sqlitedb
from zipfile import ZipFile, ZIP_DEFLATED
import time

class CAREdb(sqlitedb):
    '''
    具体的PlantCARE数据库类，封装了对数据库的基本操作
    方法：
    loadDB2dict：加载数据库到字典LOC_ID/Organism：Motif：Sequence：中
    newDB：初始化数据库（建表和索引）
    buildDB：导入指定目录下的FIMO结果文件进入数据库
    importMotifs：导入PlantCARE那里得到的Motif信息
    addEntry：Scan结果的一个条目（一行）进入数据库（查询数据库保证无重复）
    readCARE：读取PlantCARE文件中的内容至列表中
    属性：
    tables：表的结构
    views：视图的结构
    LOC_ID：LOC_ID表的内容（字典）
    Organism：Organism表的内容（字典）
    Motif：Motif表的内容（字典）
    Sequence：Sequence表的内容（字典）
    '''
    tables=[
        """
        CREATE TABLE Organism (
            id INTEGER PRIMARY KEY,
            Organism VARCHAR(50) NOT NULL UNIQUE)
        """,
        """
        CREATE TABLE Motif (
            id INTEGER PRIMARY KEY,
            Motif VARCHAR(20) NOT NULL UNIQUE,
            Description VARCHAR(120),
            Note VARCHAR(100))
        """,
        """
        CREATE TABLE Accession (
            id INTEGER PRIMARY KEY,
            Accession VARCHAR(50) NOT NULL UNIQUE)
        """,
        """
        CREATE TABLE Instance (
            id INTEGER PRIMARY KEY,
            REF_Motif INTEGER REFERENCES Motif(id),
            Type VARCHAR(12) NOT NULL,
            REF_Organism INTEGER REFERENCES Organism(id),
            Description VARCHAR(120),
            REF_Accession INTEGER REFERENCES Accession(id),
            REF_MotifSeq  INTEGER REFERENCES MotifSeq(id))
        """,
        """
        CREATE TABLE MotifSeq (
            id INTEGER PRIMARY KEY,
            Sequence VARCHAR(120) NOT NULL UNIQUE)
        """,
        """
        CREATE TABLE Scanned (
            id INTEGER PRIMARY KEY,
            REF_SeqName INTEGER REFERENCES fa_file(id),
            REF_MotifSeq INTEGER REFERENCES MotifSeq(id),
            start INTEGER NOT NULL,
            stop INTEGER NOT NULL,
            strand BOOLEAN NOT NULL,
            pValue REAL NOT NULL,
            REF_FoundMotif INTEGER REFERENCES FoundMotif(id))
        """,
        """
        CREATE TABLE fa_file (
            id INTEGER PRIMARY KEY,
            SeqName VARCHAR(20))
        """,
        """
        CREATE TABLE FoundMotif (
            id INTEGER PRIMARY KEY,
            FoundMotif VARCHAR(120) NOT NULL UNIQUE)
        """
        ]

    views=[
    """
    CREATE VIEW Main AS
    SELECT Motif,Type,Organism,Motif.Description,Accession,Sequence
    FROM Instance,Motif,Organism,Accession,MotifSeq
    WHERE Motif.id=Instance.REF_Motif
        AND Organism.id=Instance.REF_Organism
        AND Accession.id=Instance.REF_Accession
        AND MotifSeq.id=Instance.REF_MotifSeq
    """,
    """
    CREATE VIEW Motif_Seq AS
    SELECT DISTINCT Motif,Sequence
    FROM Instance,Motif,MotifSeq
    WHERE Motif.id=Instance.REF_Motif
        AND MotifSeq.id=Instance.REF_MotifSeq
    """,
    """
    CREATE VIEW ScannedMotif AS
    SELECT SeqName,Sequence AS MotifSeq,start,stop,strand,pValue,FoundMotif
    FROM Scanned,MotifSeq,fa_file,FoundMotif
    WHERE fa_file.id=Scanned.REF_SeqName
        AND MotifSeq.id=Scanned.REF_MotifSeq
        AND FoundMotif.id=Scanned.REF_FoundMotif
    """,
    """
    CREATE VIEW Motif_SeqName AS
    SELECT REF_Motif,REF_SeqName
    FROM instance,Scanned
    WHERE instance.REF_MotifSeq=Scanned.REF_MotifSeq
    """]

    Accession={}
    Organism={}
    Motif={}
    MotifSeq={}
    SeqName={}
    FoundMotif={}
    dbfile=''

    def __init__(self, dbfile):
        '''
        dbfile数据库文件路径
        '''
        self.dbfile=dbfile
        if not os.path.exists(dbfile):
            #若不存在数据库则新建之
            sqlitedb.__init__(self, dbfile)
            self.newDB(dbfile)
        else:
            sqlitedb.__init__(self, dbfile)
            self.SeqName=self.loadDB2dict('fa_file',1,0)
            self.Motif=self.loadDB2dict('Motif',1,0)
            self.MotifSeq=self.loadDB2dict('MotifSeq',1,0)
            self.FoundMotif=self.loadDB2dict('FoundMotif',1,0)
        #加载数据库四个表的内容到四个属性中
        #self.Organism=self.loadDB2dict('Organism', 1, 0)
        #self.Motif=self.loadDB2dict('Motif', 1, 0)
        #self.Accession=self.loadDB2dict('Accession', 1, 0)

    def loadDB2dict(self, tname, kcol, vcol):
        '''
        查询数据库中的表，指定两列分别作为字典的键和值
        tname：表名
        kcol：作为键的列
        vcol：作为值的列
        return 字典
        '''
        dic={}
        for rs in self.select(tname):
            dic[rs[kcol].encode('gbk')]=rs[vcol]
        return dic.copy()

    def newDB(self, dbfile):
        '''
        根据tables和views属性建立表和视图
        '''
        for SQL in self.tables+self.views:
            try:
                self.cur.execute(SQL)
            except sqlite3.OperationalError, Error:
                print "ERROR:",Error

    def importMotifs(self, care_file):
        '''
        导入PlantCARE那里得到的Motif信息
        '''
        for CARE in self.readCARE(care_file):
            if CARE[0]=='':
                #is Instance
                if CARE[4]=='':
                    #has no Sequence
                    continue
                REF_Accession=self.addEntry('Accession', \
                ['Accession'], \
                [CARE[3]], \
                self.Accession, 0)
                REF_MotifSeq=self.addEntry('MotifSeq', \
                    ['Sequence'], \
                    [multi2iupac(CARE[4])], \
                    self.MotifSeq, 0)
                self.addEntry('Instance', \
                ['REF_Motif', 'Type', 'REF_Organism', 'Description', 'REF_Accession', 'REF_MotifSeq'], \
                [REF_Motif, CARE[1], REF_Organism, CARE[2], REF_Accession, REF_MotifSeq])
            else:
                #is Motif
                REF_Organism=self.addEntry('Organism', \
                    ['Organism'], \
                    [CARE[1]], \
                    self.Organism, 0)
                REF_Motif=self.addEntry('Motif', \
                    ['Motif', 'Description', 'Note'], \
                    [CARE[0], CARE[2], CARE[3]], \
                    self.Motif, 0)
        self.commit()
        return

    def buildDBfromZip(self):
        zpath="%s/Scan" % self.dbfile.split('.')[0]
        for seq in self.MotifSeq.keys():
            print "adding %s" % seq
            scanned=self.loadScanned("%s/%s.zip" % (zpath,seq))
            self.addScanned(scanned)

    def loadScanned(self,scanned_file):
        '''
        加载扫描完的文件
        '''
        zf = ZipFile(scanned_file, 'r')
        scanned_file=zf.namelist()[0]
        scanned = [rs.split('\t') for rs in zf.read(scanned_file).split('\n')[1:-1]]
        zf.close()
        return scanned

    def addScanned(self, scan_res):
        '''
        导入FIMO的扫描结果进入数据库
        '''
        for rs in scan_res:
            (MotifSeq, SeqName, start, stop, strand, pValue, FoundMotif) = \
            (rs[0],rs[1],int(rs[2]),int(rs[3]),rs[4],float(rs[6]),rs[7])
            if strand=='+':
                strand=1
            elif strand=='-':
                strand=0
            REF_MotifSeq=self.MotifSeq[MotifSeq]
            REF_SeqName=self.addEntry('fa_file', \
                ['SeqName'], \
                [SeqName], \
                self.SeqName, 0)
            REF_FoundMotif=self.addEntry('FoundMotif', \
                ['FoundMotif'], \
                [FoundMotif], \
                self.FoundMotif, 0)
            self.addEntry('Scanned', \
                ['REF_MotifSeq', 'REF_SeqName', 'start', 'stop', 'strand', \
                'pValue', 'REF_FoundMotif'], \
                [REF_MotifSeq, REF_SeqName, start, stop, strand, \
                pValue, REF_FoundMotif])
        self.commit()
        return

    def addEntry(self, tname, colums, values, dic=None, kcol=None):
        '''
        向指定表格插入条目，用字典控制重复
        效率较高
        tname：表名
        colums：要插入的列
        values：要插入的值
        dic：要插入的表所对应的属性，如LOC_ID、Motif等
        kcol：根据values[kcol]来检查是否重复
        return 插入条目的id
        '''
        if kcol==None and dic==None:
            #即导入CARE的情况下
            self.insert(tname, colums, values, False)
        elif dic.has_key(values[kcol]):
            #有key的情况下
            return dic[values[kcol]]
        else:
            #没key就插入
            self.insert(tname, colums, values, False)
            id=len(dic)+1
            dic[values[kcol]]=id
            return id

    def readCARE(self, care_file):
        '''
        读取PlantCARE的内容至列表care中
        care_file：care文件路径
        return care（列表）
        '''
        file=open(care_file)
        care=[]
        for line in file.readlines():
            care.append(line.strip('\n').split('\t'))
        return care

    def cache(self, table='Scanned', colums=None, SeqName=None, MotifSeq=None, exact=True):
        '''
        将Scanned表读入内存，带来约三倍的速度提升。
        '''
        if SeqName or MotifSeq or (exact and table=='Scanned'):
            where={}
            if SeqName:
                where['REF_SeqName']=[str(name_id) for name_id in SeqName]
            if MotifSeq:
                where['REF_MotifSeq']=[str(motifseq_id) for motifseq_id in MotifSeq]
            if exact and table=='Scanned':
                SQL="REF_FoundMotif in (select id from FoundMotif where FoundMotif in (select Sequence from MotifSeq))"
                where['SQL']=SQL
        else:
            where=None
        self.select(table,colums,where,True)

def iupac2sites(Sequence):
    '''
    把符合IUPAC定义的多义字母序列转成具体的多条序列
    '''
    Sequences=['']
    iupac={
        'M': ['A', 'C'],
        'R': ['A', 'G'],
        'W': ['A', 'T'],
        'S': ['C', 'G'],
        'Y': ['C', 'T'],
        'K': ['G', 'T'],
        'V': ['A', 'C', 'G'],
        'H': ['A', 'C', 'T'],
        'D': ['A', 'G', 'T'],
        'B': ['G', 'C', 'T'],
        'N': ['A', 'C', 'G', 'T'],}
    for letter in Sequence:
        if iupac.has_key(letter):
            Sequences=[p+l for p in Sequences for l in iupac[letter]]
        else:
            Sequences=[p+letter for p in Sequences]
    return Sequences


def multi2sites(Sequence):
    '''
    把形如(A/C)多义字母序列转成具体的多条序列
    '''
    Sequence=Sequence.replace(")","(").split("(")
    Sequences=['']
    for letter in Sequence:
        if letter.count("/"):
            Sequences=[p+l for p in Sequences for l in letter.split("/")]
        else:
            Sequences=[p+letter for p in Sequences]
    return Sequences

def any2sites(Sequence):
    '''
    '''
    if Sequence.count('N'):
        Sequences=[]
        for seq in multi2sites(Sequence):
            Sequences.extend(iupac2sites(seq))
        return Sequences
    else:
        return multi2sites(Sequence)

def multi2iupac(Sequence):
    '''
    将形如(A/T)的序列转成IUPAC字符
    '''
    iupac={
        'M': ['A', 'C'],
        'R': ['A', 'G'],
        'W': ['A', 'T'],
        'S': ['C', 'G'],
        'Y': ['C', 'T'],
        'K': ['G', 'T']}
    Sequence=Sequence.upper()
    for key in iupac.keys():
        Sequence=Sequence.replace("("+"/".join(iupac[key])+")",key)
        Sequence=Sequence.replace("("+"/".join(iupac[key][::-1])+")",key)
    return Sequence

def writeFasta(name, Sequences):
    '''
    将序列写入sites
    '''
    file=open(name+'.sites','w')
    for seq in Sequences:
        file.write('>'+name+'\n')
        file.write(seq+'\n')
    file.close()

def writeCSV(name, data):
    file=open(name+".csv",'w')
    file.writelines(data)
    file.close()

#加载Motif
def loadMotifs(dbname):
    '''
    从数据库加载Motifs序列到字典中
    '''
    CARE=CAREdb(dbname)
    Motif_Seq={}
    for rs in CARE.select('Main'):
        Motif_name=rs[0].encode('gbk')
        Sequence=rs[1].encode('gbk')
        if not Motif_Seq.has_key(Motif_name):
            Motif_Seq[Motif_name]=[]
        if Sequence.count('/'):
            Motif_Seq[Motif_name].extend(multi2sites(Sequence))
        else:
            Motif_Seq[Motif_name].append(Sequence)
    return Motif_Seq

def loadMotifs2list(dbname):
    '''
    从数据库加载Motifs序列到列表中
    '''
    SQL="""
    SELECT DISTINCT Motif, Description, Sequence
    FROM Main
    """
    CARE=CAREdb(dbname)
    Motifs=[]
    for rs in CARE.cur.execute(SQL):
        rs=[colum.encode('gbk') for colum in rs]
        rs[-1]=multi2iupac(rs[-1])
        Motifs.append(rs)
    return Motifs

def loadScan(scanned_file):
    return [rs.split('\t') for rs in open(scanned_file).readlines()[1:-1]]

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
    dbname="%s.db" % sys.argv[1].split('.')[0]
    CARE=CAREdb(dbname)
    CARE.importMotifs('CARE.txt')
    CARE.buildDBfromZip()
    print time.time()-timestamp
