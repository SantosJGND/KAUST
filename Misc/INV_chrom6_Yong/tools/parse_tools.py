
import numpy as np 
import pandas as pd 

import collections

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)



def subset_input(genotype, summary, subset= [], bornes= 5e3):
    
    if len(subset):

        pos_sub= np.array(summary.POS,dtype= int)
        pos_sub= (pos_sub > (subset[0] - bornes)) & (pos_sub <= (subset[1] + bornes))

        summary= summary[pos_sub]
        summary= summary.reset_index(drop= True)

        genotype= genotype[:,pos_sub]
    
    return genotype, summary



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
            new_outs[nwind]= nout
            keep= []
            n_current= 0

    return new_windows, new_outs


def window_split_idx(genotype,summary,Steps= 25,window_size=100):

    window_starts= list(np.arange(0,genotype.shape[1],Steps))
    if window_starts[-1] != genotype.shape[1]:
        window_starts.append(genotype.shape[1])

    Windows= recursively_default_dict()
    Out= recursively_default_dict()

    lengths_winds= []

    for splyt in range(len(window_starts) - 1):
        IN= window_starts[splyt]
        OUT= window_starts[splyt] + window_size

        if OUT >= genotype.shape[1]:
            OUT= genotype.shape[1] - 1
        range_window= [IN,OUT]

        lengths_winds.append(OUT-IN)

        chrom= int(summary.CHROM[range_window[0]])

        start= int(summary.POS[range_window[0]])
        end= int(summary.POS[range_window[1]])

        Windows[chrom][start]= list(range(*range_window))
        Out[chrom][start]= end

        if OUT-IN < window_size:
            break

    return Windows, Out


