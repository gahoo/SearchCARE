#! /usr/bin/env python
#coding=utf-8
import os
from sqlitedb import sqlitedb

class CAREdb(sqlitedb):
    '''
    具体的PlantCARE数据库类，封装了对数据库的基本操作
    方法：
    loadDB2dict：加载数据库到字典LOC_ID/Organism：Motif：Sequence：中
    newDB：初始化数据库（建表和索引）
    buildDB：导入指定目录下的csv文件进入数据库
    addCSV：导入指定路径的csv文件进入数据库（查询数据库以保证无重复）
    addEntry：导入csv文件中的一个条目（一行）进入数据库（查询数据库保证无重复）
    addCSV2：导入指定路径的csv文件进入数据库（查询字典以保证无重复）
    addEntry：导入csv文件中的一个条目（一行）进入数据库（查询字典以保证无重复）
    readCSV：读取csv文件中的内容至列表中
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
            Sequence VARCHAR(100) NOT NULL)
        """]

    views=[
    """
    CREATE VIEW Main AS
    SELECT Motif,Type,Organism,Motif.Description,Accession,Sequence
    FROM Instance,Motif,Organism,Accession
    WHERE Motif.id=Instance.REF_Motif
    AND Organism.id=Instance.REF_Organism
    AND Accession.id=Instance.REF_Accession
    """,
    """
    CREATE VIEW Motif_Seq AS
    SELECT DISTINCT Motif,Sequence
    FROM Instance,Motif
    WHERE Motif.id=Instance.REF_Motif
    """]

    Accession={}
    Organism={}
    Motif={}

    def __init__(self, dbfile):
        '''
        dbfile数据库文件路径
        '''
        if not os.path.exists(dbfile):
            #若不存在数据库则新建之
            sqlitedb.__init__(self, dbfile)
            self.newDB(dbfile)
        else:
            sqlitedb.__init__(self, dbfile)
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

    def addCSV(self, csvfile):
        '''
        导入指定路径的csv文件
        '''
        for CARE in self.readCSV(csvfile):
            if CARE[0]=='':
                #is Instance
                if CARE[4]=='':
                    #has no Sequence
                    continue
                REF_Accession=self.addEntry('Accession', \
                ['Accession'], \
                [CARE[3]], \
                self.Accession, 0)
                self.addEntry('Instance', \
                ['REF_Motif', 'Type', 'REF_Organism', 'Description', 'REF_Accession', 'Sequence'], \
                [REF_Motif, CARE[1], REF_Organism, CARE[2], REF_Accession, CARE[4]])
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

    def readCSV(self, csvfile):
        '''
        读取csv的内容至列表csv中
        csvfile：csv文件路径
        return csv（列表）
        '''
        file=open(csvfile)
        csv=[]
        for line in file.readlines():
            csv.append(line.strip('\n').split('\t'))
        return csv

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


if __name__ == "__main__":
    #创建及加载数据库
    CARE=CAREdb('test.db')
    CARE.addCSV('CARE.txt')
    #creatSites('sites',loadMotifs('test.db'))
    for motif in loadMotifs2list('test.db'):
        print motif
    seq="(G/c)(C/a)ACCAAT(G/c)(G/c)CA(T/a)CCAAGCNNC(A/g)GAT(T/a)(T/a)N(T/g)(T/g)N(T/a)(T/a)(T/c)A"
    print multi2iupac(seq)
    #writeFasta('site',any2sites(seq))
    #print len(any2sites(seq))
