import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


def recompose_graph(G,ndes, verbose= False, show= False, fig= False,
                    dir_out= "./"):
    """
    """
    cities= {x: ndes[x,0] for x in range(ndes.shape[0])}
    #components= list(sorted(nx.connected_component_subgraphs(G), key = len, reverse=True))
    components= (G.subgraph(c) for c in nx.connected_components(G))
    components= list(components)
    complete= []
    total= []
    ix= 0

    for lt in components:

        if len(lt) == 1:
            total.append(list(lt.nodes()))
            
            continue

        ng= nx.DiGraph()
        for ed in lt.edges():
            s,t= ed
            wt= np.array(ndes[list(ed),1], dtype= int)  # edge weight as sum of contig lengths.
            wt= sum(wt)

            ng.add_edge(*ed,weigh=wt)

        nodes= ng.nodes()

        lpath= nx.dag_longest_path(ng)
        missed=[x for x in nodes if x not in lpath]
        #
        if verbose:
            if len(missed) > 0:
                
                print('{} kept. {} left out of {}'.format(len(lpath),len(missed), len(nodes)))

        complete.append(int(len(lpath) == len(nodes)))
        ###
        ###
        lpath_names= [cities[x] for x in lpath]
        nnodes= list(ng.nodes())
        H= nx.relabel_nodes(ng,cities)
        
        ##
        new_cols= decompG(H,cols_list)
        for z in list(set(new_cols)):
            
            gx= [x for x in range(len(new_cols)) if new_cols[x] == z]
            gx= [nnodes[x] for x in gx]
            total.append(gx)
        
        ###
        if fig:
            plt.figure(figsize=(9,9))
            ###
            col_map= [["blue","red"][int(x in lpath_names)] for x in H.nodes()]
            nx.draw(H,node_size = 800,width=3, node_color= new_cols, with_labels = True)
            plt.savefig(dir_out + "comp{}_n{}_nsub{}.png".format(ix,len(lt),len(set(new_cols))))
            #
            if show:
                plt.show()
            plt.close()

        ix+=1
    
    return complete, total 



def te_stats(edges,TEgff, TEcomp=80):
    tagged= []

    for edix in range(edges.shape[0]):
        assess= []
        for edu in [1,2]:
            ctg= edges.iloc[edix]["C{}".format(edu)]
            s0=edges.iloc[edix]["s{}".format(edu)]
            e0= np.sum(np.array([edges.iloc[edix]["s{}".format(edu)], edges.iloc[edix]["Fl{}".format(edu)]], dtype= int ))

            overlaps= TEgff[(TEgff.contig == ctg) & (TEgff.start<e0) & (TEgff.end > s0)]
            rL= edges.iloc[edix]["rL".format(edu)]
            sumOverlap=0
            prop= 0
            many= 0
            TeAv= 0
            if overlaps.shape[0] > 0:
                #print(edges.loc[edix].shape)
                #print(overlaps)

                sumSs= overlaps.start - s0
                sumSs= np.array(sumSs)
                sumEs= np.array(overlaps.end - s0)
                sumEs[sumEs > e0] = e0
                sumSs[sumSs < 0]= 0            
                sumOverlap= sumEs - sumSs
                tels= np.array(overlaps.end - overlaps.start)
                teps= 100 * sumOverlap / tels
                #print(teps)
                many= np.sum(teps > TEcomp)
                TeAv= np.median(teps)
                sumOverlap= np.sum(sumOverlap)
                prop= 100 * sumOverlap / rL



            assess.extend([sumOverlap,overlaps.shape[0], prop, many, TeAv])

        tagged.append(assess)

    tagged= pd.DataFrame(tagged)
    tagged.columns=  ["Ol1","On1","Op1","Ot1","Om1","Ol2","On2","Op2","Ot2","Om2"]

    # Ol = total overlap with contigs
    # On = Number of contigs overlapped ; 
    # Op = Proportion of fragment covered by TEs ; 
    # Ot = total number of Complete TEs (> TEcomp %)
    # Om = Average percentage of coverage across TEs 
    
    return tagged



def read_edges(edgefl, nderef, filterpd= [],filter_dict= {}, minL= 3):
    
    edge_dict= {}
    
    if len(filterpd) > 0:
        for ft,val in filter_dict.items():
            print(edgefl.shape)
            #print(filterpd.shape)
            for ent in [1,2]:
                #print(filterpd[ft.format(ent)] < val)
                edgefl= edgefl[(filterpd[ft.format(ent)] < val)].reset_index(drop=True)
                filterpd= filterpd[(filterpd[ft.format(ent)] < val)].reset_index(drop=True)
    
    edgel= np.array(edgefl, dtype= str)
    print(edgel.shape)
    print("##")
    for ix in range(edgel.shape[0]):
        inn= edgel[ix][1]
        innix= nderef.index(inn)
        out= edgel[ix][2]
        outix= nderef.index(out)

        if innix in edge_dict.keys():
            if outix not in edge_dict[innix]:
                edge_dict[innix].append(outix)
        else:
            edge_dict[innix]= [outix]

    edge_dict= {x:g for x,g in edge_dict.items() if len(g) <= minL}
    
    return edge_dict


def decompG(H, cols_list):
    Hnodes= list(H.nodes())
    Hprox= list(Hnodes)
    Vamp=H.copy()
    new_cols= ["blue"] * len(Hnodes)
    d=0
    while Hprox:
        Lpath= nx.dag_longest_path(Vamp)

        for td in Lpath:
            new_cols[Hnodes.index(td)]= str(d)#cols_list[d]
            Vamp.remove_node(td)
        nnodes= Vamp.nodes()
        singles= list(nx.isolates(Vamp))

        for td in singles:
            Vamp.remove_node(td)
        Hprox= [x for x in nnodes if x not in singles]
        
        d+=1
    
    return new_cols



############
############

import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--dir', type= str, default= "./")
parser.add_argument('--edgel', type= str, default= "reads.edges")
parser.add_argument('--nodes', type= str, default= "ref_contiglens.tsv")
parser.add_argument('--minL', type= int, default= 10)
parser.add_argument('--gff', type= str, default= "EDTA.complete.gff")
parser.add_argument('--Op', type= int, default= 40)
parser.add_argument('--Ot', type= int, default= 2)
parser.add_argument('--Om', type= int, default= 70)

args = parser.parse_args()

cols_list= ["red","purple","green","yellow","hotpink", "brown", "magenta","cyan","orange"]

############

sofdir= args.dir
mdir= sofdir
edge_filelist= args.edgel.split(",")
edge_filelist= [x for x in edge_filelist if len(x) > 0]
minL=args.minL
nodef= args.nodes

with open(nodef,"r") as fp:
    ndes= fp.readlines()
    ndes= [x.strip().split() for x in ndes]
    ndes= [[x[0],x[-1]] for x in ndes]

ndes= np.array(ndes,dtype= str)

###################################
###################################
#edges1f, edges2f= edge_filelist

edge_series= [pd.read_csv(x, sep=" ", header= None) for x in edge_filelist]

#edges2= pd.read_csv(edges2f, sep=" ", header= None)

edges= pd.concat(tuple(edge_series),ignore_index= True)
edges.columns= ["R","C1","C2","rL","c1L","c2L","Q1","Q2","s1","s2","Fl1","Fl2"]

edges["end1"]=edges.s1+edges.Fl1
edges["end2"]=edges.s2+edges.Fl1

print(edges.shape)

#edge_dict= read_edges(edge_filelist, list(ndes[:,0]), minL= minL)

###

TEedta= args.gff
TEgff= pd.read_csv(TEedta,sep= "\t",header= None)
TEgff.columns= ["contig","source","feature","start","end","score","strand","frame","att"]
TEgff.head()

TEcomp=80

tagged= te_stats(edges,TEgff, TEcomp=TEcomp)
tagged.to_csv(mdir + "edgeTE_stats.tsv")
#####################################
filter_dict= {
    "Op{}": args.Op,
    "Ot{}": args.Ot,
    "Om{}": args.Om
}

minL= args.minL

print(tagged.head())
edge_dict= read_edges(edges, list(ndes[:,0]), tagged, filter_dict, minL= minL)


############################################

G= nx.Graph()

my_nodes=range(ndes.shape[0])
G.add_nodes_from(my_nodes)

for s,clan in edge_dict.items():
    for e in clan:
        G.add_edge(s,e)

#####
######

dir_out= mdir + "comps_plots/"
from pathlib import Path
Path(dir_out).mkdir(parents=True, exist_ok=True)

complete, total= recompose_graph(G,ndes, verbose= False, show= False,
                    dir_out= "./")


###### ## Calculate summary_stat
sizes= [[int(ndes[x,1]) for x in y] for y in total]
sizes= sorted([sum(x) for x in sizes], reverse= True)
total_size= sum(sizes)

L50= 0
so= 0
while so < total_size / 2:
    so += sizes[L50]
    L50 +=1

N50= sizes[L50]

summ= pd.DataFrame([[len(total),L50, N50]], columns= ["totalC","L50","N50"])

summ.to_csv(mdir + "summary.txt", index= False, sep= "\t")

