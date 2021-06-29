import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles

def filter_info(invbed, accessions, dist_threshold= 5e3,query="query"):
    bed_filtered= []
    print(invbed.chrom.unique())
    for chrom in invbed.chrom.unique():
        bedsel= invbed.loc[invbed.chrom==chrom].reset_index(drop= True)
        bedsel["IN"]= pd.to_numeric(bedsel['IN'])
        bedsel["OUT"]= pd.to_numeric(bedsel['OUT'])
        bedsel= bedsel.sort_values(by="IN").reset_index(drop= True)
        proxy_bed= bedsel.copy()
        #print(proxy_bed.head(50))

        INV_labels= ins_labels(proxy_bed, dist_threshold= dist_threshold)

        proxy_bed["INV_merge"]= merge_labs(INV_labels)
        #print(proxy_bed.head(50))
        bedsel= bed_collapse(proxy_bed,accessions,query=query)
        bedsel["IN"]= pd.to_numeric(bedsel["IN"])
        bedsel["OUT"]= pd.to_numeric(bedsel["OUT"])
        print(bedsel.head(30))
        print("#######")
        bed_filtered.append(bedsel)

    bed_filtered= pd.concat(bed_filtered)
    return bed_filtered


from sklearn.cluster import AgglomerativeClustering
def ins_labels(bedsel,inv_tab= ["IN","OUT"],dist_threshold= 1000.0):
    """

    """
    print(bedsel.head())
    inv_arrays= [np.array(bedsel[y],dtype= int) for y in inv_tab]
    inv_labels= []
    print(bedsel.shape)
    for idx in range(len(inv_tab)):
        inv_distances= inv_arrays[idx]
        
        clustering = AgglomerativeClustering(n_clusters= None,distance_threshold= dist_threshold)

        clusters= clustering.fit(inv_distances.reshape(-1,1))
        inv_labels.append(clusters.labels_)

    ovlab=[0]
    dfm= bedsel.iloc[0]["OUT"]
    d=0
    ix= 1
    print( bedsel.iloc[ix]["IN"],bedsel.iloc[ix-1]["OUT"], dist_threshold)
    print("h")
    print( bedsel.iloc[ix-1]["OUT"])
    for ix in range(1,bedsel.shape[0]):
        lt= bedsel.iloc[ix]["OUT"]
        if (bedsel.iloc[ix]["IN"] - dist_threshold) < dfm:
            lt= max([lt, dfm])
            #
        else: d+=1
        
        #if bedsel.iloc[ix]["IN"] > (bedsel.iloc[ix-1]["OUT"] + dist_threshold):
        #    d+=1
        dfm= lt
        ovlab.append(d)

    inv_labels.append(np.array(ovlab))

    inv_labels= np.array(inv_labels).T
    print(inv_labels.shape)
    print(inv_labels[:5])

    return inv_labels


def merge_labs(label_array):
    '''
    merge array of label columns.
    '''

    unique_now= label_array[0,:]
    nlabs= []
    labi= 0

    for idx in range(label_array.shape[0]):
        comp_now= label_array[idx,:]
        comparison= comp_now == unique_now
        compar= sum((comparison[0],comparison[1]))

        if compar < 2 and comparison[2] == False:
            labi += 1
        
	
        nlabs.append(labi)
        

        unique_now= comp_now

    return nlabs



def bed_collapse(bed_array,accessions,clcol= 'INV_merge',start_cols= ['chrom','IN'],
                 end_cols= ['OUT'],func_do= [np.min, np.max],query="query"):
    '''
    collapse array based on labels column.
    '''
    #print(bed_array.head)
    inv_clusts= list(bed_array[clcol])
    inv_dict= {z: [x for x in range(bed_array.shape[0]) if inv_clusts[x] == z] for z in set(inv_clusts)}
    clust_set= set(inv_clusts)
    new_array= []
    new_cols= []
    new_cols.extend(start_cols)
    new_cols.extend(end_cols)

    new_cols.extend(['L'] + accessions)
    #new_cols.append('idx')
    for cl,g in inv_dict.items():

        new_row= []
        minmax= []
        acc_dict= {x: [] for x in accessions}

        for col in start_cols:
            ends= [int(bed_array.iloc[x][col]) for x in g]
            
            if col == 'IN':
                ends= int(func_do[0](ends))
                minmax.append(ends)
            else:
                ends= ends[0]

            new_row.append(ends)

        for col in end_cols:
            ends= [int(bed_array.iloc[x][col]) for x in g]

            if col == 'OUT':
                ends= int(func_do[1](ends))
                minmax.append(ends)
            else:
                ends= ends[0]

            new_row.append(ends)

        if len(minmax)== 2:
            le= minmax[1] - minmax[0]

            new_row.append(int(le))

        clumped=[]

        for clp in g:
            acc= bed_array.iloc[clp][query]
            acc_dict[acc].append(bed_array.iloc[clp]["SV"].strip(" "))

        for acc in accessions:
            if len(acc_dict[acc]) == 0:
                clumped.append("NA")
            else:
                clumped.append(",".join(acc_dict[acc]))
        new_row.extend(clumped)
        #clumped= ";".join(clumped)
        #print(bed_array["index1"])
        #new_row.append(".".join([str(bed_array.iloc[x]["index1"]) for x in g]))
        #new_row.append(".".join([str(x) for x in g]))
        #new_row.append(';'.join([str(x) for x in clumped]))
        new_array.append(new_row)
    #print(clumped)
    new_array= np.array(new_array)
    #print(new_array[:4])
    #print(new_cols)
    new_array= pd.DataFrame(new_array,columns= new_cols)
    return new_array


#######################################################################
######################################################################

import sys
import itertools as it

def import_split(flist,snames,id='s1',pavs= {},minl= 50, maxl= 1e6):
        outdict= {x: [] for x in pavs.keys()}
        for ix in range(len(flist)):
                fl=flist[ix]
                sn= snames[ix]
                #print(ix)
                #print(fl)
                with open(fl,'r') as fp:
                        ln= fp.readline()
                        
                        while ln:
                                #ln= ln.replace('Chr','')
                                ln= ln.replace('Chr0','')
                                ln= ln.replace('Chr','')

                                for pt,tels in pavs.items():
                                        if len([x for x in tels if x in ln]) > 0 and "NC_" not in ln and "NW_" not in ln:
                                                outdict[pt].append([ln.strip().split()[x] for x in [0,1,2,4,3]] + [id,sn])
                                ln= fp.readline()

        for ix in outdict.keys():
                outdict[ix]= pd.DataFrame(outdict[ix])
                print(ix)                
                print(outdict[ix].shape)
                outdict[ix].columns= ["chrom", "IN", "OUT", "L", "SV", "acc", "soft"]
                outdict[ix][["IN", "OUT","L"]].apply(pd.to_numeric)
                outdict[ix]["L"]= np.array(outdict[ix]["L"], dtype= int)
                outdict[ix]= outdict[ix][(outdict[ix]["L"] > minl) & (outdict[ix]["L"] < maxl)]
                outdict[ix]= outdict[ix].drop(columns=["L"])
        return outdict


def plot_bars(pavarray, softs,percent= False, fig= "", show= True):
    pav_stats= []

    pav_pres= np.sum(pavarray,axis= 1)

    pav_pres= np.sum(pavarray,axis= 1)
    print(pavarray[:50])
    soft_total= np.sum(pavarray,axis= 0)
    unique, counts = np.unique(pav_pres, return_counts=True)
    perc= counts[0] / sum(counts)
    pav_stats.append(soft_total)

    if percent:
        counts= [100 * x / sum(counts) for x in counts]
    pav_stats.append(counts)

    ###
    pavrep= pavarray[pav_pres == 1] # print(pavarray[pav_pres == 1][:10])
    counts= np.sum(pavrep,axis= 0) # np.argmax(pavarray[pav_pres == 1],axis= 1)
    print(pavrep[:10],softs, pavarray[:5])
    #pavrep= [softs[x] for x in pavrep]
    #unique, counts = np.unique(pavrep, return_counts=True)
    #resh= [softs.index(x) for x in unique]
    #counts= np.array([counts[x] for x in resh])
    print("total")
    print(counts, soft_total)
    if percent:
        soft_unique= 100 * counts / soft_total
        counts= [100 * x / np.sum(counts) for x in counts]

    pav_stats.extend([counts, soft_unique])

    ###
    titles= ["soft_total","sfs","unique","soft_unique"]
    labs= ["soft total","sfs (freq= {})".format(percent),"unique (freq= {})".format(percent),"soft_unique"]
    for ix in range(len(pav_stats)):
        plt.figure()
        ir= pav_stats[ix]
        rug= softs
        if titles[ix] == "sfs":
            rug= range(1,len(rug)+1)
        
        plt.bar(rug, pav_stats[ix])
        plt.title(titles[ix])
        plt.xticks(fontsize= 13)
        plt.yticks(fontsize= 13)
        plt.xlabel("",fontsize= 13)
        plt.ylabel(labs[ix],fontsize= 13)
        if fig:
            plt.savefig(fig + titles[ix] + "_bar.png")
        if show:
            plt.show()

        plt.close()
    
    pav_stats= pd.DataFrame(np.array(pav_stats).T, columns= titles)
    pav_stats["soft"]= softs
    return pav_stats


def plot_bars1(pavarray, softs,percent= False, fig= "", show= True):
    
    pav_pres= np.sum(pavarray,axis= 1)
    print(pav_pres[:10])
    soft_total= np.sum(pavarray,axis= 0)
    unique, counts = np.unique(pav_pres, return_counts=True)
    perc= counts[0] / sum(counts)
    print(soft_total)
    
    plt.figure()
    plt.ylabel("Count", fontsize= 13)
    
    if percent:
        counts= [100 * x / sum(counts) for x in counts]
        
        plt.ylabel("%", fontsize= 13)

    plt.title('{} % unique'.format(perc))
    plt.xticks(fontsize= 13)
    plt.yticks(fontsize= 13)
    plt.xlabel("",fontsize= 13)
    plt.bar(unique, counts)

    if fig:
        plt.savefig(fig + "Dist_bar.png")
    if show:
        plt.show()

    ########
    ########
    pavrep= np.argmax(pavarray[pav_pres == 1],axis= 1)
    print(len(pavrep), pavarray.shape)
    pavrep= [softs[x] for x in pavrep]
    unique, counts = np.unique(pavrep, return_counts=True)

    plt.figure()
    if percent:
        counts= counts / soft_total
        print(counts)
        plt.ylabel("%", fontsize= 13)

    plt.bar(unique, counts)
    plt.title("unique dist.")
    plt.xticks(fontsize= 13)
    plt.yticks(fontsize= 13)
    plt.xlabel("",fontsize= 13)
    plt.ylabel("",fontsize= 13)
    if fig:
        plt.savefig(fig + "Soft_PercUniqe_bar.png")
    if show:
        plt.show()
    plt.close()

    unique, counts = np.unique(pavrep, return_counts=True)

    plt.figure()
    if percent:
        counts= [100 * x / sum(counts) for x in counts]
        plt.ylabel("%", fontsize= 13)

    plt.bar(unique, counts)
    plt.title("unique dist.")
    plt.xticks(fontsize= 13)
    plt.yticks(fontsize= 13)
    plt.xlabel("",fontsize= 13)
    plt.ylabel("",fontsize= 13)
    if fig:
        plt.savefig(fig + "unique_bar.png")
    if show:
        plt.show()

    plt.close()
    


def stats_comp(pavar, fig='', show= '',cols= []):
    pav_pres= np.sum(pavar,axis= 1)
    if len(cols) == 0:
        cols= ['s{}'.format(x) for x in range(pavar.shape[1])]

    for comb in it.combinations(list(range(pavar.shape[1])), 3):
        npav= pavar[:,list(comb)]
        combs= []
        
        nnames= [cols[x] for x in comb]
        print(nnames)
        nset= []
        for ix in range(3):
            nset.append(npav[(np.sum(npav,axis=1) == 1) & (npav[:,ix] == 1),:].shape[0])

        for cmp in list(it.combinations([0,1,2], 2)):
            print(cmp)
            nnames.append('-'.join([cols[x] for x in cmp]))

            nset.append(npav[(np.sum(npav,axis=1) == 2) & (np.sum(npav[:,cmp],axis=1)== 2),:].shape[0])

        nset.append(npav[np.sum(npav,axis=1) == 3].shape[0])
        print(nset)
        print(np.sum(nset[:2]))
        nnames.append("all")


        plt.figure(figsize=(10,10))
        nset=np.array(nset)[[0,1,3,2,4,5,6]]
        venn3(subsets=tuple(nset), set_labels= nnames)
        if fig:
            plt.savefig(fig + "threeway_venn_{}.png".format("-".join([cols[x] for x in comb])))
        if show:
            plt.show()

        plt.close()


