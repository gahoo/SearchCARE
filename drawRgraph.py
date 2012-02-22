#! /usr/bin/env python
#coding=utf-8
import rpy2.robjects as robj
from rpy2.robjects.packages import importr
from rpy2.rinterface import RRuntimeError
import sys
import os
import time
from zipfile import ZipFile, ZIP_DEFLATED
import cPickle

def drawRHist(Motif_Pos,outpath):
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    xrange=robj.IntVector([-3000,0])
    for motif in Motif_Pos.keys():
        print motif
        if not Motif_Pos[motif]:
            print "empty"
            continue
        dist=-robj.IntVector(Motif_Pos[motif]).ro
        motif=motif.replace('/','_')
        grdevices.png(file="%s/%s.png" % (outpath,motif), width=512, height=512)
        graphics.hist(dist,breaks=50,main=motif,xlab="Distribution",xlim=xrange)
        #robj.r('hist(0-%s, xlab="x", main="hist(x)")' %dist.r_repr())
        grdevices.dev_off()

def drawRDensity(Motif_Pos,outpath):
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    xrange=robj.IntVector([-3000,0])
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
        graphics.plot(density,main=motif,xlab="Distribution",xlim=xrange)
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
    zf = ZipFile('%s.zpkl' % basename, 'r')
    Motif_Pos = cPickle.loads(zf.open('motif_dist.pkl').read())
    zf.close()
    drawRHist(Motif_Pos,outpath="%s/hist" % basename)
    drawRDensity(Motif_Pos,outpath="%s/density" % basename)
    drawRmultidensity(Motif_Pos,outpath=basename,top=5)
    print time.time()-timestamp
