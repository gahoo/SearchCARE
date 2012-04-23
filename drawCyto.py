#! /usr/bin/env python
#coding=utf-8
import sys
import os
import time
import getopt
import cPickle
import xmlrpclib

'''
读取Stats导出的Node和Edge，通过CytoscapeRPC导出Cytoscape图片和gml格式网络
'''

def loadEdges(edgefile):
    '''
    加载Edge文件
    返回字典Edges
    '''
    edges=open(edgefile,'r')
    (motifAs,motifBs,co_nums,A_nums,B_nums,Ps)=([],[],[],[],[],[])
    for edge in edges.readlines():
        (motifA,motifB,co_num,A_num,B_num,p)=edge.split('\t')
        motifAs.append(motifA)
        motifBs.append(motifB)
        co_nums.append(int(co_num))
        A_nums.append(int(A_num))
        B_nums.append(int(B_num))
        Ps.append(float(p))
    Edges={'motifA':motifAs,'motifB':motifBs,'co_nums':co_nums,'A_nums':A_nums,'B_nums':B_nums,'Ps':Ps}
    return Edges

def loadNodes(nodefile):
    '''
    加载Node文件
    返回字典Nodes
    '''
    nodes=open(nodefile,'r')
    (motifs,descs,lst_cnts,db_cnts,lst_sizes,db_sizes,enrich_ratios,Ps)=([],[],[],[],[],[],[],[])
    for node in nodes.readlines():
        (motif,desc,lst_cnt,db_cnt,lst_size,db_size,enrich_ratio,p)=node.split('\t')
        motifs.append(motif)
        descs.append(desc)
        lst_cnts.append(int(lst_cnt))
        db_cnts.append(int(db_cnt))
        lst_sizes.append(int(lst_size))
        db_sizes.append(int(db_size))
        enrich_ratios.append(float(enrich_ratio))
        Ps.append(float(p))
    Nodes={'motifs':motifs,'descs':descs,'lst_cnts':lst_cnts,'db_cnts':db_cnts,'lst_sizes':lst_sizes,'db_sizes':db_sizes,'enrich_ratios':enrich_ratios,'Ps':Ps}
    return Nodes

def createNetwork(server,Nodes,Edges,network_name,layout,enrichratioThresh,pThresh):
    '''
    根据Node和Edge构建网络
    '''
    networkid = server.Cytoscape.createNetwork(network_name)
    importNodes(server,networkid,Nodes)
    importEdges(server,networkid,Edges)
    server.Cytoscape.setDefaultBackgroundColor('default',"#FFFFFF")
    server.Cytoscape.performLayout(networkid,layout)
    setMappers(server,Nodes,Edges)
    saveNet(server,networkid,network_name)
    if enrichratioThresh:
        subNetEnrichNodes(server,networkid,network_name,Nodes,layout,ER=enrichratioThresh)
    if pThresh:
        subNetSigEdges(server,networkid,network_name,Nodes,Edges,layout,pThresh)
    server.Cytoscape.destroyNetwork(networkid)

def importNodes(server,networkid,Nodes):
    '''
    导入节点属性信息
    '''
    server.Cytoscape.addStringNodeAttributes('descs', Nodes['motifs'], Nodes['descs'])
    server.Cytoscape.addIntegerNodeAttributes('lst_cnts', Nodes['motifs'], Nodes['lst_cnts'])
    server.Cytoscape.addIntegerNodeAttributes('db_cnts', Nodes['motifs'], Nodes['db_cnts'])
    server.Cytoscape.addIntegerNodeAttributes('lst_sizes', Nodes['motifs'], Nodes['lst_sizes'])
    server.Cytoscape.addIntegerNodeAttributes('db_sizes', Nodes['motifs'], Nodes['db_sizes'])
    server.Cytoscape.addDoubleNodeAttributes('enrich_ratios', Nodes['motifs'], Nodes['enrich_ratios'])
    server.Cytoscape.addDoubleNodeAttributes('Ps', Nodes['motifs'], Nodes['Ps'])
    #节点类型，是否富集，用于设置节点形状
    nodes_types=map(lambda x: 'enrich' if x<=0.05 else 'not', Nodes['Ps'])
    server.Cytoscape.addStringNodeAttributes('nodes_types', Nodes['motifs'], nodes_types)

def importEdges(server,networkid,Edges):
    '''
    导入边属性信息
    '''
    server.Cytoscape.createEdges(Edges['motifA'],Edges['motifB'])
    edgenames=["%s (directed) %s" % (Edges['motifA'][i],Edges['motifB'][i]) for i in range(len(Edges['motifB']))]
    server.Cytoscape.addIntegerEdgeAttributes('co_nums',edgenames,Edges['co_nums'])
    server.Cytoscape.addIntegerEdgeAttributes('A_nums',edgenames,Edges['A_nums'])
    server.Cytoscape.addIntegerEdgeAttributes('B_nums',edgenames,Edges['B_nums'])
    server.Cytoscape.addDoubleEdgeAttributes('Ps',edgenames,Edges['Ps'])
    edges_types=map(lambda x: 'sign' if x<=0.01 else 'not', Edges['Ps'])
    server.Cytoscape.addStringEdgeAttributes('edges_types', edgenames, edges_types)

def setMappers(server,Nodes,Edges):
    '''
    将属性映射至图像中，包括颜色，大小，形状等
    节点：三角形：富集、大小：出现次数、颜色：富集率
    边：实线：显著相关P<0.01、宽度：P显著程度、颜色：共现次数
    '''
    #Nodes Attributes
    #富集的节点为三角形，否则为原型
    server.Cytoscape.createDiscreteMapper('default','nodes_types','Node Shape', \
                                          "ellipse", \
                                          {'enrich':"triangle",'not':"ellipse"})
    #节点大小与lst_cnts，即列表中模体出现次数成正比
    server.Cytoscape.createContinuousMapper('default','lst_cnts','Node Size', \
                                            [float(min(Nodes['lst_cnts'])), float(max(Nodes['lst_cnts']))], \
                                            [5.0,10.0,50.0,60.0])
    #节点颜色深浅与enrich_ratios，即富集率成正比
    server.Cytoscape.createContinuousMapper('default','enrich_ratios', 'Node Color', \
                                            [float(min(Nodes['enrich_ratios'])), 1.0, float(max(Nodes['enrich_ratios']))], \
                                            ['#FFFFFF', '#FFFFFF', '#FFFFFF', '#0000FF', '#0000FF'])
    #Edges Attributes
    #显著相关的边为实线，一般相关为虚线
    server.Cytoscape.createDiscreteMapper('default','edges_types','Edge Line Style', \
                                          "EQUAL_DASH", \
                                          {'sign':"SOLID",'not':"EQUAL_DASH"})
    #边的宽度与P值成反比，P越小线越粗
    server.Cytoscape.createContinuousMapper('default','Ps','Edge Line Width', \
                                            [float(min(Edges['Ps'])), float(max(Edges['Ps']))], \
                                            [8.0,7.0,0.5,0.3])
    #边的颜色和共现次数co_nums成正比，颜色越深共现次数越多
    server.Cytoscape.createContinuousMapper('default','co_nums', \
                                            'Edge Color', \
                                            [float(min(Edges['co_nums'])), float(max(Edges['co_nums']))], \
                                            ['#FAFAFF', '#FAFAFF', '#0000FF', '#0000FF'])

def saveNet(server,networkid,network_name):
    '''
    保存网络到文件，导出为图片
    '''
    #计算合适的scale
    node_cnt=server.Cytoscape.countNodes(networkid)
    scale=2+node_cnt/15.0
    #避免motif名称出边界
    zoom=0.8*server.Cytoscape.getZoom(networkid)
    server.Cytoscape.setZoom(networkid,zoom)
    netfile="%s/%s_Cyto" % (os.getcwd(),network_name)
    print network_name,scale,zoom#,netfile
    #server.Cytoscape.saveNetwork(networkid,"%s.gml" % netfile)
    server.Cytoscape.executeCommand('network','export',{'file':"%s.xgmml" % netfile,'type':"xgmml"})
    server.Cytoscape.exportView(networkid,"%s.png" % netfile,"png",scale)

def subNetEnrichNodes(server,networkid,network_name,Nodes,layout,ER=1.0):
    '''
    节点富集率大于指定值的子网络
    '''
    enriched=[Nodes['motifs'][i] for i in range(len(Nodes['enrich_ratios'])) if Nodes['enrich_ratios'][i]>=ER]
    #去掉那些虽然enrich但不在网络中的节点
    all_nodes=server.Cytoscape.getNodes(networkid)
    enriched=list(set(enriched)&set(all_nodes))
    server.Cytoscape.selectNodes(networkid,enriched)
    sub_networkid=server.Cytoscape.createNetworkFromSelection(networkid,"%s Nodes ER>%s" % (network_name,str(ER)))
    server.Cytoscape.performLayout(sub_networkid,layout)
    if server.Cytoscape.countNodes(sub_networkid):
        saveNet(server,sub_networkid,"%s_ER%s" % (network_name,str(ER)))
    server.Cytoscape.destroyNetwork(sub_networkid)

def subNetSigEdges(server,networkid,network_name,Nodes,Edges,layout,pThresh):
    '''
    边p值小于指定值的子网络
    '''
    #由于不能直接选择边创建网络，只能迂回
    #先根据显著的边所连的所有节点创建一个临时网络
    significant_edges=["%s (directed) %s" % (Edges['motifA'][i],Edges['motifB'][i]) for i in range(len(Edges['Ps'])) if Edges['Ps'][i]<=pThresh]
    server.Cytoscape.selectEdges(networkid,significant_edges)
    selectNodesBYEdges(server,networkid,significant_edges)
    tmp_networkid=server.Cytoscape.createNetworkFromSelection(networkid,"tmp")
    #再删除临时网络中不显著的边
    sub_edges=server.Cytoscape.getEdges(tmp_networkid)
    sub_edges_p=server.Cytoscape.getEdgesAttributes('Ps',sub_edges)
    not_sign_edges=[sub_edges[i] for i in range(len(sub_edges_p)) if sub_edges_p[i]>pThresh]
    #not_sign_edges=map(lambda x: True if x<=0.001 else False, sub_edges_p)
    for not_sign_edge in not_sign_edges:
        server.Cytoscape.removeEdge(tmp_networkid,not_sign_edge)
    #最后选择剩余所有边所连接的节点，据此创建网络
    remaining_edges=server.Cytoscape.getEdges(tmp_networkid)
    selectNodesBYEdges(server,tmp_networkid,remaining_edges)
    sub_networkid=server.Cytoscape.createNetworkFromSelection(tmp_networkid,"%s Edges P<%s" % (network_name,str(pThresh)))
    server.Cytoscape.performLayout(sub_networkid,layout)
    if server.Cytoscape.countNodes(sub_networkid):
        saveNet(server,sub_networkid,"%s_EdgeP%s" % (network_name,str(pThresh)))
    server.Cytoscape.destroyNetwork(tmp_networkid)
    server.Cytoscape.destroyNetwork(sub_networkid)

def selectNodesBYEdges(server,networkid,edges):
    '''
    根据所选的边选择对应的节点
    '''
    targer_nodes=server.Cytoscape.getEdgeTargetNodes(networkid,edges)
    source_nodes=server.Cytoscape.getEdgeSourceNodes(networkid,edges)
    server.Cytoscape.selectNodes(networkid,targer_nodes)
    server.Cytoscape.selectNodes(networkid,source_nodes)

def usage(server):
    print "drawCyto.py -e <edgefile> -n <nodefile> [-m] network_name [-l] layout [-r] enrichratioThresh [-p] pThresh"
    #print server.Cytoscape.getNodeShapeNames()
    #print server.Cytoscape.getLineStyleNames()
    print "Layouts: %s" % ",".join(server.Cytoscape.getLayoutNames())

def has_file(filename):
    if not os.path.exists(filename):
        raise Exception,"There is no %s here" % filename

if __name__ == '__main__':
    server = xmlrpclib.ServerProxy("http://localhost:9000")
    #init var
    (network_name, layout, enrichratioThresh, pThresh)=(None,"kamada-kawai", 1.0, 0.001)
    try:
        server.Cytoscape.test()
    except:
        print "CytoscapeRPC is not running"
        sys.exit(2)
    try:
        opts, args = getopt.getopt(sys.argv[1:], "he:n:m:l:r:p:", ["help","edgefile=","nodefile=","network_name=","layout=","enrichratioThresh=","pThresh="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(server)
            sys.exit()
        elif opt in ("-e", "--edgefile"):
            edgefile=arg
            has_file(edgefile)
            if not network_name:
                network_name=edgefile.split('.')[0]
        elif opt in ("-n", "--nodefile"):
            nodefile=arg
            has_file(nodefile)
            if not network_name:
                network_name=nodefile.split('.')[0]
        elif opt in ("-m", "--network_name"):
            network_name=arg
        elif opt in ("-l", "--layout"):
            layout=arg
        elif opt in ("-r", "--enrichratioThresh"):
            enrichratioThresh=float(arg)
        elif opt in ("-p", "--pThresh"):
            pThresh=float(arg)
    if not 'edgefile' or not 'nodefile' in dir():
        print "tell me both edgefile and nodefile please"
        usage(server)
        sys.exit(2)

    print network_name,layout,enrichratioThresh,pThresh
    timestamp=time.time()
    createNetwork(server, \
                loadNodes(nodefile), \
                loadEdges(edgefile), \
                network_name, \
                layout, \
                enrichratioThresh, \
                pThresh)
    print time.time()-timestamp
