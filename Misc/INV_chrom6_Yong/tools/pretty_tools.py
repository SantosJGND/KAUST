import numpy as np
from structure_tools.vcf_geno_tools import read_geno_nanum
import pandas as pd

def proc_name(Names):
    '''locally specific de-dedoubling names.'''
    for x in range(len(Names)):
        ind= Names[x]
        newid= ind.split('_')
        if len(newid) > 2:
            newid= '_'.join(newid[:2])
        else:
            newid= newid[0]

        Names[x]= newid
    
    return Names


def pretty_input(Home= 'D:/GitHub/Tools_and_toys/VCF_analysis/Extract/vcf/', info_file= '3K_info.txt',Chr= 6, row_info= 6, header_info= 9,
                                                                                phased= False):

    filename= Home + 'Extract_Chr' + str(Chr) + '_15000.vcf'

    row_info= 6
    header_info= 9
    phased= False

    genotype, summary, Names= read_geno_nanum(filename, row_info= row_info, header_info= header_info,phased= phased)

    ## Process Names vcf names.
    ## Instance specific processing due to ID copy in VCF file.
    name_proc= 1
    if name_proc:
        Names= proc_name(Names)
        name_proc= 0


    print('Number of markers: {}'.format(genotype.shape[1]))
    print('Number of individuals: {}'.format(genotype.shape[0]))

    RG_info= pd.read_csv(info_file,sep= '\t')
    
    return genotype, summary, Names, RG_info


################################################
################################################
from sklearn.decomposition import PCA
import plotly

import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

from tools.vcf_geno_tools import geno_subset_random

def plot_interest(gen_sample, subsummary, RG_info, subset_col, Names, code= {},ref_keep= [0,1,2],
                                                    gp_names_dict= {}, ID_col= 'IRIS_ID',
                                                    include= [], Sn= 400, Sm= 5000, gp_cols= {}, sample= True):
    
    if sample:
        gen_sample, subsummary, code_vec, code_lib, Nsample, Msample= geno_subset_random(gen_sample,summary, RG_info, 
                                                                                         ID_col,subset_col, Names,code=code, include= Names_select, 
                                                                                         Sn= Sn, Sm= Sm)
    else:
        Nsample= list(range(len(Names)))

    PCsel= 1

    Names_subsel= [Names[x] for x in Nsample]
    ###


    gp_cols= {z:'rgb({})'.format(','.join([str(x) for x in g])) for z,g in gp_cols.items()}

    ref_dict= list(RG_info['IRIS_ID'])
    ref_dict= [ref_dict.index(x) for x in Names_subsel]
    ref_code= [RG_info['Initial_subpop'][x] for x in ref_dict]

    ref_dict= [['admx',x][int(x in code.keys())] for x in ref_code]
    ref_dict= [code[x] for x in ref_dict]

    keep_vec= [x for x in range(len(ref_dict)) if ref_dict[x] in ref_keep or Names_subsel[x] in Names_select]
    names_kept= [Names_subsel[x] for x in keep_vec]
    Names_idx= [names_kept.index(x) for x in Names_select]
    ref_dict= [ref_dict[x] for x in keep_vec]


    ref_dict= {z: [x for x in range(len(ref_dict)) if ref_dict[x]==z] for z in set(ref_dict)}
    ###
    pca = PCA(n_components=3)

    pca.fit(gen_sample[keep_vec])

    pc_data= pca.transform(gen_sample[keep_vec])

    ###

    fig= [go.Scatter3d(
        x= pc_data[g,0],
        y= pc_data[g,PCsel],
        z= pc_data[g,PCsel+1],
        text= [Names_subsel[x] for x in g],
        mode= "markers",
        name= gp_names_dict[z],
        marker= dict(
            size= 5,
            color= gp_cols[z]
        )
    ) for z,g in ref_dict.items()]

    fig.append(go.Scatter3d(
        x= pc_data[Names_idx[:5],0],
        y= pc_data[Names_idx[:5],PCsel],
        z= pc_data[Names_idx[:5],PCsel+1],
        mode= "markers",
        name= 'interest subtrop',
        text= Names_select[:5],
        marker= dict(
            size= 5,
            color= 'rgb(30,144,255)',
            symbol= 'x',
            line= dict(width= 1)
        )
        )
    )

    fig.append(go.Scatter3d(
        x= pc_data[Names_idx[5:10],0],
        y= pc_data[Names_idx[5:10],PCsel],
        z= pc_data[Names_idx[5:10],PCsel+1],
        mode= "markers",
        name= 'control subtrop',
        text= Names_select[5:10],
        marker= dict(
            size= 5,
            color= 'rgb(255,215,0)',
            symbol= 'cross',
            line= dict(width= 2)
        )
        )
    )

    fig.append(go.Scatter3d(
        x= pc_data[Names_idx[10:],0],
        y= pc_data[Names_idx[10:],PCsel],
        z= pc_data[Names_idx[10:],PCsel+1],
        mode= "markers",
        name= 'control temp',
        text= Names_select[10:],
        marker= dict(
            size= 5,
            color= 'rgb(173,255,47)',
            symbol= 'diamond',
            line= dict(width= 2)
        )
    ))



    layout= go.Layout()

    Figure= go.Figure(data= fig,layout= layout)
    iplot(Figure)

