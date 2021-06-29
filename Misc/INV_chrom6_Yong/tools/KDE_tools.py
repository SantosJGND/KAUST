
import numpy as np
import pandas as pd
import re
import scipy
import itertools as it

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import MeanShift, estimate_bandwidth

import collections

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


#### reference KDE

from tools.vcf_tools import samp_sample, samp_same, samp_same_v2c1, samp_same_v2c2


def ms_clustersget(data,Bandwidth_split= 20):

    Accurate = []
    params = {'bandwidth': np.linspace(np.min(data), np.max(data),Bandwidth_split)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0,cv=3)
    
    ######################################
    ####### TEST global Likelihood #######
    ######################################
    Focus_labels = range(data.shape[0])
    
    #### Mean Shift approach
    ## from sklearn.cluster import MeanShift, estimate_bandwidth
    
    bandwidth = estimate_bandwidth(data[Focus_labels,:], quantile=0.2, n_samples=len(Focus_labels))
    if bandwidth <= 1e-3:
        bandwidth = 0.1
    
    ms = MeanShift(bandwidth=bandwidth, cluster_all=False, min_bin_freq=25)
    ms.fit(data[Focus_labels,:])
    labels = ms.labels_
    
    Tree = {x:[Focus_labels[y] for y in range(len(labels)) if labels[y] == x] for x in [g for g in list(set(labels)) if g != -1]}

    SpaceX = {x:data[Tree[x],:] for x in Tree.keys()}
    
    
    Construct= {}

    for hill in SpaceX.keys():
        if len(Tree[hill]) <= 3:
            continue
        grid.fit(data[Tree[hill],:])
        
        # use the best estimator to compute the kernel density estimate
        kde = grid.best_estimator_
        
        P_dist = kde.score_samples(data[Tree[hill],:])
        Dist = kde.score_samples(data)
        P_dist= np.nan_to_num(P_dist)
        Dist= np.nan_to_num(Dist)
        if np.std(P_dist) == 0:
            Dist= np.array([int(Dist[x] in P_dist) for x in range(len(Dist))])
        else:
            Dist = scipy.stats.norm(np.mean(P_dist),np.std(P_dist)).cdf(Dist)
        Dist= np.nan_to_num(Dist)
        Construct[hill] = Dist

    return Construct
        
    


def extract_profiles(global_data,genotype,seq_idx, ref_labels,
                    RG_info, ID_col,subset_col,Names, ref_lib= {}, n_comps= 4,
                    repn= 100, code= {},others= 'admx',Sn= 500, same= True,clean= False,
                    ms_comp= False):
    '''
    Extract KDE profiles for specific accessions (global_idx) from reference groups in PCA space. 
    Reduction and KDE calculated at seq_idx positions in genotype array.
    Reference accessions from ref_labels groups are permuted. samp_sample() function is used to 
    sample accessions using RG_info, takes Sm
    '''
    ## estimate the bandwith

    pca2 = PCA(n_components=n_comps, whiten=False,svd_solver='randomized')

    cluster_profiles= {x:[] for x in ref_labels}
    
    ## perform KDE.
    combine= {}
    tkeys= ref_labels
    var_comp_store= []

    construct= {}

    for rp in range(repn):
        print(rp)
        if same:
            Names_idx, kde_class_labels, kde_label_dict, Nsample= samp_same_v2c2(ref_lib,Names,Sn= Sn)
        else:
            Names_idx, kde_class_labels, kde_label_dict, Nsample= samp_sample(genotype,RG_info, ID_col,subset_col,Names,code= code,others= others,Sn= Sn)

        dat_foc= genotype[:,seq_idx]
        dat_foc= dat_foc[global_data]

        Sequences= genotype[:,seq_idx]
        Sequences= Sequences[Names_idx]

        if Sequences.shape[1] <= 3:
            Results[Chr][c] = [0,0]
            print('hi')
            continue

        pca2.fit(Sequences)
        data= pca2.transform(Sequences)
        data_ref= pca2.transform(dat_foc)

        local_pcvar= list(pca2.explained_variance_ratio_)

        #local_pcvar= [local_pcvar]

        var_comp_store.append(local_pcvar)

        params = {'bandwidth': np.linspace(np.min(data), np.max(data),15)}
        grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,cv= 3,verbose=0)
        
        ref_q= []

        for bull in tkeys:
            
            Quanted_set= data[kde_label_dict[bull],:]
            grid.fit(Quanted_set)
            kde = grid.best_estimator_

            P_dist = kde.score_samples(Quanted_set)
            Fist = kde.score_samples(data_ref)

            ## Normalizing log-likelihood estimates by those of the reference set and extracting their cdf.
            #Fist = scipy.stats.norm(np.mean(P_dist),np.std(P_dist)).cdf(Fist)

            cluster_profiles[bull].append(Fist)

        ####
        ####
        construct= ms_clustersget(data)


    cluster_profiles= {
        x: np.array(g) for x,g in cluster_profiles.items()
    }

    cluster_profiles= {
        x: np.median(g,axis= 0) for x,g in cluster_profiles.items()
    }

    var_comp_store= np.array(var_comp_store)
    var_comp_store= np.median(var_comp_store,axis= 0)

    return cluster_profiles, var_comp_store, construct


def KDE_window_profiles(genotype,Windows,target_idx,
    RG_info, ID_col,subset_col,Names,ref_labels,Chr= 1,n_comps= 4, exclude= [], gp_focus= [],
    repn= 100, code= {},others= 'admx',Sn= 500, Sm= 10000, same= True, clean= False,ms_comp= False):
    '''
    extract reference profiles across windows. 
    '''

    Windows_profiles= recursively_default_dict()
    construct_lib= {Chr:{}}
    var_comp_store= []

    ref_lib= {}
    if same:
        ref_lib= samp_same_v2c1(genotype,RG_info, ID_col,subset_col,Names,exclude= exclude,code= code,others= others)

    if gp_focus:
        ref_lib= {z:g for z,g in ref_lib.items() if z in gp_focus}

    for bl in Windows[Chr].keys():

        print("window {}".format(bl))

        seq_idx= Windows[Chr][bl]

        profiles, var_comps, construct= extract_profiles(target_idx,genotype,seq_idx,ref_labels,
                                                RG_info, ID_col,subset_col,Names, ref_lib= ref_lib, n_comps= n_comps,
                                                repn= repn, code= code,others= others,Sn= Sn, same= same,clean= clean,
                                                ms_comp= ms_comp)

        ### store stuff.

        Windows_profiles[Chr][bl]= profiles
        construct_lib[Chr][bl]= construct
        var_comps= [bl, *list(var_comps)]
        var_comp_store.append(var_comps)


    var_comp_store= np.array(var_comp_store)
    var_comp_store= pd.DataFrame(var_comp_store,columns=['set',*['PC' + str(x + 1) for x in range(n_comps)]])
    if ms_comp:
        return Windows_profiles, var_comp_store, construct_lib
    else:
        return Windows_profiles, var_comp_store

