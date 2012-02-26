#! /usr/bin/env python
#coding=utf-8
import rpy2.robjects as robj
from rpy2.robjects.packages import importr
from rpy2.rinterface import RRuntimeError
import sys
import os
import time
import getopt
from zipfile import ZipFile, ZIP_DEFLATED
import cPickle

def drawRHist(Motif_Pos,outpath):
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    xlim=robj.IntVector([-3000,0])
    for motif in Motif_Pos.keys():
        print motif
        if not Motif_Pos[motif]:
            print "empty"
            continue
        dist=-robj.IntVector(Motif_Pos[motif]).ro
        motif=motif.replace('/','_')
        grdevices.png(file="%s/%s.png" % (outpath,motif), width=512, height=512)
        graphics.hist(dist,breaks=50,main=motif,xlab="Distribution",xlim=xlim)
        #robj.r('hist(0-%s, xlab="x", main="hist(x)")' %dist.r_repr())
        grdevices.dev_off()

def drawRDensity(Motif_Pos,outpath):
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    xlim=robj.IntVector([-3000,0])
    for motif in Motif_Pos.keys():
        print motif
        if not Motif_Pos[motif]:
            print "empty"
            continue
        dist=-robj.IntVector(Motif_Pos[motif]).ro
        try:
            density=robj.r.density(dist)
        except RRuntimeError:
            continue
        motif=motif.replace('/','_')
        grdevices.png(file="%s/%s.png" % (outpath,motif), width=512, height=512)
        graphics.plot(density,main=motif,xlab="Distribution",xlim=xlim)
        grdevices.dev_off()

def drawRmultidensity(Motif_Pos,outpath,col=False,top=None):
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    geneplotter = importr('geneplotter')
    RColorBrewer = importr('RColorBrewer')
    grdevices.png(file="%s.png" % outpath, width=512, height=512)
    Pos={}
    if not top:
        top=len(Motif_Pos)
    for motif in Motif_Pos.keys():
        Pos[motif]=-robj.IntVector(Motif_Pos[motif]).ro
    Pos=robj.ListVector(Pos)
    names=robj.r('names(sort(sapply(%s,length),decreasing=T))' % Pos.r_repr())
    if col:
        colors=robj.r('colorRampPalette(brewer.pal(9,"BrBG"))')
        colors=robj.r.rev(colors(top))
        geneplotter.multidensity(Pos.rx(names[:top]),col=colors,lwd=3,xlab="Location")
    else:
        geneplotter.multidensity(Pos.rx(names[:top]),lwd=3,xlab="Location")
    grdevices.dev_off()

def drawRGraph(motif_dist,outpath):
    '''
    motif_dist[motif][MotifSeq]={'+':dist,'-':dist}
    '''
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    geneplotter = importr('geneplotter')
    RColorBrewer = importr('RColorBrewer')
    xlim=robj.IntVector([-3000,0])
    ylim=robj.FloatVector([0,0.002])
    for motif in motif_dist.keys():
        print motif
        if not motif_dist[motif]:
            print "empty"
            continue
        motif_seq={}
        #画multidensity
        for MotifSeq in motif_dist[motif]:
            #merge+,-
            pos=mergeList(motif_dist[motif][MotifSeq])
            if len(pos)>1:
                motif_seq[MotifSeq]=-robj.IntVector(pos).ro
        #merge MotifSeq
        dist=mergeList(motif_seq)
        #print motif_seq
        if not motif_seq:
            print "empty"
            continue
        motif_seq=robj.ListVector(motif_seq)
        names=robj.r('names(sort(sapply(%s,length),decreasing=T))' % motif_seq.r_repr())
        #print motif_seq.r_repr(),names.r_repr()
        motif=motif.replace('/','_')
        grdevices.png(file="%s/%s.png" % (outpath,motif), width=512, height=512)
        dist=robj.IntVector(dist)
        density=robj.r.density(dist)
        if len(motif_seq)>1:
            geneplotter.multidensity(motif_seq.rx(names),lwd=3,xlab="Distribution",main=motif,xlim=xlim,ylim=ylim)
        else:
            graphics.plot(density,lty='dashed',lwd=3,main=motif,xlab="Distribution",xlim=xlim,ylim=ylim)
        #graphics.hist(dist,add=True,breaks=50)
        graphics.rug(dist)
        graphics.lines(density,lty='dashed',lwd='4')
        grdevices.dev_off()
    return

def mergeList(motif_dict):
    '''
    motif_dict={'...':[],'...':[]}
    '''
    pos=[]
    [pos.extend(v) for v in motif_dict.values()]
    return pos

def drawRbarplot(SeqName_counts,output):
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    #SeqName_counts=robj.ListVector(SeqName_counts)
    counts=robj.IntVector(SeqName_counts.values())
    counts.names=robj.StrVector(SeqName_counts.keys())
    #还有排序方面的问题
    grdevices.png(file="%s/BarPlot.png" % output, width=512, height=512)
    graphics.barplot(counts)
    grdevices.dev_off()


def usage():
    print "CAREdb.py -f <zpkl> [-o] <output>"

def has_file(filename):
    if not os.path.exists(filename):
        raise Exception,"There is no %s here" % filename

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:o:", ["help","filename=","output="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-f", "--filename"):
            filename=arg
            has_file(filename)
            output=filename.split('.')[0]
        elif opt in ("-o", "--output"):
            output=arg
    if not 'filename' in dir():
        print "tell me filename please"
        usage()
        sys.exit(2)
    if not os.path.exists(output):
        os.makedirs(output)
    timestamp=time.time()
    zf = ZipFile(filename, 'r')
    motif_dist = cPickle.loads(zf.open('motif_dist.pkl').read())
    SeqName_counts = cPickle.loads(zf.open('SeqName_counts.pkl').read())
    zf.close()
    #drawRHist(Motif_Pos,outpath="%s/hist" % basename)
    #drawRDensity(Motif_Pos,outpath="%s/density" % basename)
    #drawRmultidensity(Motif_Pos,outpath=basename,top=5)
    #drawRGraph(motif_dist,output)
    drawRbarplot(SeqName_counts,output)
    print time.time()-timestamp