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

def drawRGraph(motif_dist,enriched,outpath):
    '''
    motif_dist[motif][MotifSeq]={'+':dist,'-':dist}
    '''
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    geneplotter = importr('geneplotter')
    RColorBrewer = importr('RColorBrewer')
    xlim=robj.IntVector([-3000,0])
    ylim=robj.FloatVector([0,0.002])
    for motif in enriched.keys():
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

def drawRbarplot(SeqName_counts,enriched,output,top=25):
    if not enriched:
        return
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    #SeqName_counts=robj.ListVector(SeqName_counts)
    counts=robj.IntVector(SeqName_counts.values())
    counts.names=robj.StrVector(SeqName_counts.keys())
    #counts=robj.r.sort(counts)
    p=robj.FloatVector([enriched[motif]['p'] for motif in enriched.keys()])
    p.names=robj.StrVector(enriched.keys())
    #enriched_names=robj.StrVector(enriched.keys())
    p=robj.r.sort(p,decreasing=True)
    enriched_counts=counts.rx(p.names)
    #top_counts=robj.r.tail(counts,n=top)
    #print top_counts.r_repr()
    grdevices.png(file="%s/%s_bar.png" % (output,output), width=512, height=512)
    margin=robj.IntVector([3,9,4,2])
    graphics.par(mar=margin)
    bar=graphics.barplot(enriched_counts,main="Enriched motifs counts",horiz=True,las=1,col='lightblue')
    graphics.text(x=enriched_counts,y=bar,label=robj.r.signif(p,digits=2),po=2) 
    #graphics.text(bar,labels=top_counts,pos=4,offset=10)
    grdevices.dev_off()

def drawRboxplot(Merged_dist,enriched,output):
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    #Merged_dist=dict([(k,list(v)) for k,v in Merged_dist.items()])
    for motif in Merged_dist.keys():
        Merged_dist[motif]=-robj.IntVector(Merged_dist[motif]).ro
    Merged_dist=robj.ListVector(Merged_dist)
    names=robj.r('names(sort(sapply(%s,length),decreasing=T))' % Merged_dist.r_repr())
    enriched_names=robj.StrVector(enriched.keys())
    grdevices.png(file="%s/%s_box.png" % (output,output), width=512, height=512)
    margin=robj.IntVector([3,9,4,2])
    graphics.par(mar=margin)
    graphics.boxplot(Merged_dist.rx(enriched_names),main="Boxplot of Motif Positions",horizontal=True,las=1,col='lightblue')
    grdevices.dev_off()

def drawRCoocur(co_Occur,outpath):
    '''
    co_Occur[(motifa,motifb)][REF_SeqName]={(motifa_pos,motifb_pos):distance}
    '''
    for motif_pair in co_Occur.keys():
        print motif_pair,len(co_Occur[motif_pair])
        drawRdistance(co_Occur[motif_pair],outpath,motif_pair)
        drawRCoDistibution(co_Occur[motif_pair],outpath,motif_pair)

def drawRdistance(motif_pair,outpath,pair_name):
    distances=[]
    filename="&".join(pair_name).replace('/','_')
    for REF_SeqName in motif_pair.keys():
        distances.extend(motif_pair[REF_SeqName].values())
    dist=robj.IntVector(distances)
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    grdevices.png(file="%s/%s_distance.png" % (outpath,filename), width=512, height=512)
    graphics.hist(dist,breaks=50,border=4,main="Distance Histogram of \n%s" % "&".join(pair_name),xlab="Distances")
    grdevices.dev_off()

def drawRCoDistibution(motif_pair,outpath,pair_name):
    motifA_pos=[]
    motifB_pos=[]
    filename="&".join(pair_name).replace('/','_')
    (motifa,motifb)=pair_name
    for REF_SeqName in motif_pair.keys():
        #print motif_pair[REF_SeqName].keys()
        for motifa_pos,motifb_pos in motif_pair[REF_SeqName].keys():
            motifA_pos.append(motifa_pos)
            motifB_pos.append(motifb_pos)
            #print motifa_pos,motifb_pos
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    geneplotter = importr('geneplotter')
    motifA_pos=-robj.IntVector(motifA_pos).ro
    motifB_pos=-robj.IntVector(motifB_pos).ro
    Pos={motifa:motifA_pos,motifb:motifB_pos}
    Pos=robj.ListVector(Pos)
    grdevices.png(file="%s/%s_distribution.png" % (outpath,filename), width=512, height=512)
    geneplotter.multidensity(Pos.rx(),lwd=3,xlab="Distribution",main="Distribution of \n%s" % filename)
    graphics.rug(motifA_pos,col=4)
    graphics.rug(motifB_pos,col=2)
    grdevices.dev_off()
    #Scatter plot
    grdevices.png(file="%s/%s_Scatter.png" % (outpath,filename), width=512, height=512)
    limit=robj.IntVector([-3000,0])
    graphics.plot(motifA_pos,motifB_pos,main="Position Scatter Plot of\n%s&%s" % (motifa,motifb), \
                    xlab="Positions of %s" % motifa, \
                    ylab="Positions of %s" % motifb, \
                    xlim=limit,ylim=limit)
    graphics.abline(1,1)
    grdevices.dev_off()



def usage():
    print "CAREdb.py -z <zpkl> [-o] <output>"

def has_file(filename):
    if not os.path.exists(filename):
        raise Exception,"There is no %s here" % filename

if __name__ == '__main__':
    co_dist_range=(10,300)
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hz:o:t:r:", ["help","zpklfile=","output=","top=","co_dist_range="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-z", "--zpklfile"):
            zpklfile=arg
            has_file(zpklfile)
            output=zpklfile.split('.')[0]
            top=25
        elif opt in ("-o", "--output"):
            output=arg
        elif opt in ("-t", "--top="):
            top=arg
        elif opt in ("-r", "--co_dist_range"):
            #arg=lower#higer
            co_dist_range=tuple(map(int,arg.split("#")))
    if not 'zpklfile' in dir():
        print "tell me filename please"
        usage()
        sys.exit(2)
    if not os.path.exists(output):
        os.makedirs(output)
    if not os.path.exists("%s/Enriched" % output):
        os.makedirs("%s/Enriched" % output)
    if not os.path.exists("%s/CoOccur" % output):
        os.makedirs("%s/CoOccur" % output)
    timestamp=time.time()
    zf = ZipFile(zpklfile, 'r')
    Merged_MotifSeq_dist = cPickle.loads(zf.open('Merged_MotifSeq_dist.pkl').read())
    SeqName_counts = cPickle.loads(zf.open('SeqName_counts.pkl').read())
    Merged_dist=cPickle.loads(zf.open('Merged_dist.pkl').read())
    enriched=cPickle.loads(zf.open('Enriched.pkl').read())
    try:
        co_Occur=cPickle.loads(zf.open('co_Occur_%s#%s.pkl' % co_dist_range).read())
    except KeyError:
        co_Occur=None
    zf.close()
    #drawRHist(Motif_Pos,outpath="%s/hist" % basename)
    #drawRDensity(Motif_Pos,outpath="%s/density" % basename)
    #drawRmultidensity(Motif_Pos,outpath=basename,top=5)
    drawRGraph(Merged_MotifSeq_dist,enriched,"%s/Enriched" % output)
    drawRbarplot(SeqName_counts,enriched,output,top)
    drawRboxplot(Merged_dist,enriched,output)
    if co_Occur:
        drawRCoocur(co_Occur,"%s/CoOccur" % output)
    print time.time()-timestamp
