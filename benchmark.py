#! /usr/bin/env python
#coding=utf-8
import sqlite3
import sys
from sqlitedb import sqlitedb
from CAREdb import CAREdb
import time

SQL=["attach ':memory:' as cache",
    "create table cache.Scanned as select * from Scanned",
    "create table cache.Instance as select * from Instance",
    #"select REF_SeqName from Scanned where REF_MotifSeq=498",
    "select REF_SeqName from Scanned where REF_MotifSeq=498",
    """
    select distinct count(*) from
    (select REF_SeqName from Scanned where REF_MotifSeq=498
    INTERSECT
    select REF_SeqName from Scanned where REF_MotifSeq=237)
    """,
    """
    select count(distinct REF_SeqName) from Scanned where REF_MotifSeq in (
     select REF_MotifSeq from instance where REF_Motif=498)
    """,
    """
    select count(distinct REF_SeqName) from Scanned where REF_MotifSeq in (
     select distinct REF_MotifSeq from instance where REF_Motif=498)
    """,
    "select REF_MotifSeq from instance where REF_Motif=498",
    "select distinct REF_MotifSeq from instance where REF_Motif=498",
    "select count(REF_MotifSeq) from instance where REF_Motif=498",
    "select count(distinct REF_MotifSeq) from instance where REF_Motif=498",
    """
            SELECT DISTINCT REF_SeqName
            FROM cache.Scanned
            WHERE REF_MotifSeq IN (
                    SELECT REF_MotifSeq
                    FROM cache.Instance
                    WHERE REF_Motif=3)
            """,
    """
            SELECT DISTINCT REF_SeqName
            FROM Scanned
            WHERE REF_MotifSeq IN (
                    SELECT REF_MotifSeq
                    FROM Instance
                    WHERE REF_Motif=3)
            """

    ]

SQL2=["attach ':memory:' as cache",
    "create table cache.Scanned as select * from Scanned",
    "create table cache.Instance as select * from Instance",
    #"select REF_SeqName from cache.Scanned where REF_MotifSeq=498",
    """
    select distinct count(*) from
    (select REF_SeqName from Scanned where REF_MotifSeq=1
    INTERSECT
    select REF_SeqName from Scanned where REF_MotifSeq=3)
    """,
    """
    select count(distinct REF_SeqName) from
    (select REF_SeqName from cache.Scanned where REF_MotifSeq=1
    INTERSECT
    select REF_SeqName from cache.Scanned where REF_MotifSeq=3)
    """,
    """
    select count(distinct REF_SeqName) from
    (select distinct REF_SeqName from cache.Scanned where REF_MotifSeq in (
     select REF_MotifSeq from cache.Instance where REF_Motif=1)
    INTERSECT
    select distinct REF_SeqName from cache.Scanned where REF_MotifSeq in (
         select REF_MotifSeq from cache.Instance where REF_Motif=3))
    """,
    """
    select count(distinct REF_SeqName) from
    (select distinct REF_SeqName from cache.Scanned where REF_MotifSeq in (
     select REF_MotifSeq from Instance where REF_Motif=1)
    INTERSECT
    select distinct REF_SeqName from cache.Scanned where REF_MotifSeq in (
         select REF_MotifSeq from Instance where REF_Motif=3))
    """,
    """
    select count(distinct REF_SeqName) from
    (select distinct REF_SeqName from cache.Scanned where REF_MotifSeq in (
     select distinct REF_MotifSeq from cache.Instance where REF_Motif=1)
    INTERSECT
    select distinct REF_SeqName from cache.Scanned where REF_MotifSeq in (
         select distinct REF_MotifSeq from cache.Instance where REF_Motif=3))
    """,
    """
    select count(distinct REF_SeqName) from cache.Scanned where REF_MotifSeq in (
     select REF_MotifSeq from cache.Instance where REF_Motif=1)
    """,
    """
    select count(distinct REF_SeqName) from cache.Scanned where REF_MotifSeq in (
     select REF_MotifSeq from cache.Instance where REF_Motif=1)
    """,
    """
    select count(distinct REF_SeqName) from cache.Instance,cache.Scanned
    where REF_Motif=1 and cache.Scanned.REF_MotifSeq=cache.Instance.REF_MotifSeq
    """,
    """
    select count(distinct REF_SeqName) from cache.Instance,cache.Scanned
    where cache.Scanned.REF_MotifSeq=cache.Instance.REF_MotifSeq and REF_Motif=1
    """
    ]

SQL3=[
    """
    SELECT REF_SeqName,REF_Motif,REF_Organism,Description,start,stop,strand,pValue
    FROM cache.Instance,cache.Scanned
    WHERE cache.Instance.REF_MotifSeq=cache.Scanned.REF_MotifSeq
    AND REF_SeqName in ('40747','10669','1620','5828','30256','10402')
    """,
    ]

benchmark=CAREdb('RAP_3kbp.db')
benchmark.cache('Instance')
benchmark.cache('Scanned',['REF_SeqName','REF_MotifSeq','start','stop','strand','pValue'],SeqName=['40747','10669','1620','5828','30256','10402'])

for sql in SQL3:
    print sql
    timestamp=time.time()
    benchmark.cur.execute(sql)
    print [rs[0] for rs in benchmark.cur.fetchall()]
    print time.time()-timestamp
'''
print "------------"
tables=['Scanned','cache.Scanned']
for tname in tables:
    timestamp=time.time()
    #benchmark.select(tname)
    benchmark.cache(colums=['REF_SeqName','REF_MotifSeq'])
    print time.time()-timestamp

benchmark=CAREdb('RAP_3kbp.db')

timestamp=time.time()
#benchmark.select('Motif_SeqName',['REF_SeqName','REF_Motif'])
print time.time()-timestamp
print "--------------"
timestamp=time.time()
#benchmark.cache('Scanned',['REF_SeqName','REF_MotifSeq','start'])
#benchmark.select('cache.Motif_SeqName',['REF_SeqName','REF_Motif'])
print len(benchmark.loadDB2dict('fa_file',1,0))
print time.time()-timestamp
'''