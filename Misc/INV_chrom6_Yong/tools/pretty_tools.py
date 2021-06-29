import numpy as np
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


from tools.input_utilities import read_vcf_allel



def pretty_input(vcf_file= '', info_file= '3K_info.txt', haps_extract= False,ploidy= 2, ID_col= "IRIS_ID"):
    
    if not len(vcf_file): 
        print('no vcf file provided.')

        return ''
    genotype, summary, Names= read_vcf_allel(vcf_file,haps_extract= haps_extract)   


    ## Process Names vcf names.
    ## Instance specific processing due to ID copy in VCF file.
    
    Names= proc_name(Names)

    Nnames= []

    if haps_extract:
        #Names= list(Names) + list(Names)
        for idx in range(ploidy):
            new_names= [x + '_{}'.format(idx) for x in Names]
            Nnames.extend(new_names)

        Names= list(Nnames)
    else:
        Names= list(Names)

    print('Number of markers: {}'.format(genotype.shape[1]))
    print('Number of individuals: {}'.format(genotype.shape[0]))

    RG_info= pd.read_csv(info_file,sep= '\t')

    RG_store= []
    rgNames_store= []
    
    if haps_extract:
        for idx in range(ploidy):
            RG_store.append(RG_info)
            rgNames_store.extend([x + "_{}".format(idx) for x in list(RG_info[ID_col])])

        RG_info= pd.concat(RG_store)
        RG_info[ID_col]= rgNames_store
        RG_info= RG_info.reset_index()
    
    return genotype, summary, Names, RG_info



