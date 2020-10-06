
import numpy as np 
import pandas as pd 


def geno_Lwind_split(summary, geno_bornes= [],Steps= 25e3,window_size=5e4):
    '''
    split genotype array into windows by length, steps.
    assumes genotype has a single chrom.
    '''
    
    POS= np.array(summary.POS,dtype= int)
    if not geno_bornes:
        geno_bornes= [0,max(POS)]
    
    window_starts= np.array(np.arange(geno_bornes[0],geno_bornes[1],Steps),dtype=int)
    
    Windows= {}
    Out= {z: z+window_size for z in window_starts}
    
    current_winds= []
    
    for idx in range(len(POS)):
        posh= int(POS[idx])
        
        if len(window_starts):
            if window_starts[0]<= posh:
                current_winds.append(window_starts[0])
                Windows[window_starts[0]]= []
                window_starts= window_starts[1:]
        
        current_rm= []
        for windx in current_winds:
            if posh > Out[windx]:
                current_rm.append(windx)
                continue
            
            Windows[windx].append(idx)
        
        current_winds= [x for x in current_winds if x not in current_rm]
            
    return Windows, Out



def wind_compress(Windows,Out,min_snp= 3):
    """
    merge contiguous windows to reach minimum number of snps.
    """

    wind_sort= sorted(Windows.keys())

    new_windows= {}
    new_outs= {}

    keep= []
    n_current= 0

    for idx in range(len(wind_sort)):
        wind= wind_sort[idx]

        keep.append(wind)
        n_current+= len(Windows[wind])

        if n_current >= min_snp:
            nwind= keep[0]
            nout= Out[keep[-1]]
            poss= [Windows[x] for x in keep]
            poss= list(it.chain(*poss))

            new_windows[nwind]= poss
            new_outs= nout
            keep= []
            n_current= 0

    return new_windows, new_outs
