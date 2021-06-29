import numpy as np 
import os
import itertools as it
from scipy.stats import kstest,poisson 
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


from scipy.stats import gaussian_kde

def draw_comp(n= 50, rangepos= [0,1], rug= 200,
              comp= 'uniform',nrep= 1000):
    
    stat_store= []
    
    
    for idx in range(nrep):
        pos= np.random.randint(*rangepos,size= n)
        
        xs = np.linspace(chrom_start,chrom_end,rug)
        density = gaussian_kde(pos)
        
        Y= density(xs)
        
        kt= kstest((pos - rangepos[0]) / rangepos[1], comp)
        
        stat_store.append([
            *kt,
            density.factor,
            np.mean(Y),
            np.std(Y)
        ])
    
    
    return np.array(stat_store)



def draw_comp_ks(n= 50, rangepos= [0,1], rug= 200,
              comp= 'uniform',nrep= 1000,
                geno_bornes= [0,1000],Steps= 5e5,window_size=1e6):
    
    stat_store= []
    count_store= []
    
    for idx in range(nrep):
        pos= np.random.randint(*rangepos,size= n)
        
        xs = np.linspace(chrom_start,chrom_end,rug)
        density = gaussian_kde(pos)
        
        Y= density(xs)
        
        kt= kstest((pos - rangepos[0]) / rangepos[1], comp)
        
        stat_store.append([
            *kt,
            density.factor,
            np.mean(Y),
            np.std(Y)
        ])
        
        ##
        ##
        pos= np.random.randint(*geno_bornes,size= n)
        pos= sorted(pos)
        #
        td= np.array(pos,dtype= int)
        td= pd.DataFrame(td,columns= ["POS"])
        
        Windows,Out= geno_Lwind_split(td, geno_bornes= [chrom_start,chrom_end],Steps= ksteps,window_size=kwindow)
        #        
        
        kp= [len(set(x)) for x in Windows.values()]

        count_store.extend(kp)
    
    
    return np.array(stat_store), count_store



def geno_Lwind_split(summary, geno_bornes= [],Steps= 25e3,window_size=5e4):
    '''
    split genotype array into windows by length, steps.
    assumes genotype has a single chrom.
    '''
    
    POS= np.array(summary.POS,dtype= int)
    if not geno_bornes:
        ## assume that chromosome does not end at last INV position; 
        ## without genome sizes, best bet is that INVs are uniformily distributed.
        geno_bornes= [0,(max(POS) + min(POS))]
    
    window_starts= np.array(np.arange(geno_bornes[0],geno_bornes[1],Steps),dtype=int)
    
    Windows= {x: [] for x in window_starts}
    Out= {z: z+window_size for z in window_starts}
    
    current_winds= []
    
    for idx in range(len(POS)):
        posh= int(POS[idx])
        
        if len(window_starts):
            if window_starts[0]<= posh:
                current_winds.append(window_starts[0])
                
                window_starts= window_starts[1:]
                
            d= 1-int(len(window_starts) > 0)
            ids= 0
            while d == 0:
                if Out[window_starts[ids]]< posh:
                    ids += 1
                else:
                    current_winds.append(window_starts[0])
                    
                    window_starts= window_starts[1:]
                    d+= 1
            
            window_starts= window_starts[ids:]
        
        current_rm= []
        
        for windx in current_winds:
            if posh > Out[windx]:
                if idx == len(POS) - 1:
                    current_rm.append(windx)
                else:
                    if POS[idx+1] > Out[windx]:
                        current_rm.append(windx)
                        
                continue
            
            Windows[windx].append(idx)
        
        current_winds= [x for x in current_winds if x not in current_rm]
            
    return Windows, Out






import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--input', type= str,
                    default="filtered_HQ_INV.bed")


parser.add_argument('--kwindow', type=int,
                    default=2e5)

parser.add_argument('--ksteps', type=int,
                    default=1e5)

parser.add_argument('--alpha', type=float,
                    default=1e-2)

parser.add_argument('--cmet', type=str,
                    default="fdr_bh")

parser.add_argument('--centro', type=str,
                    default="")


args = parser.parse_args()



kwindow= args.kwindow
ksteps= args.ksteps
alpha= args.alpha
corr_method= args.cmet


file= args.input

###
if args.centro:
    centromeres= pd.read_csv(args.centro,header= None,sep= '\t')
    centromeres.columns= ['ID','chrom','IN','OUT','L']

###
invbed= pd.read_csv(file, header=None, sep= '\t')
invbed.columns= ['chrom','IN','OUT','L']

## replace names in bed column
#invbed.chrom= invbed.chrom.str.split('.',n= 1, expand= True)[1]

invbed.head()

out_dir= file.split('/')[-1].split('.')[0] + '_{}_{}'.format(alpha,corr_method) + "/"
print(out_dir)

os.mkdir(out_dir)

#######################
#######################
pdout= []
kt_out= []
bed_filtered= []
for chrom in invbed.chrom.unique():

    bedsel= invbed.loc[invbed.chrom==chrom].reset_index(drop= True)
    bedsel= bedsel.sort_values(by="IN")

    bedsel["IN"]= pd.to_numeric(bedsel["IN"])
    bedsel["OUT"]= pd.to_numeric(bedsel["OUT"])

    print(bedsel.head(30))
    bed_filtered.append(bedsel)

    chrom_start= 0
    chrom_end= bedsel.OUT.max() + bedsel.IN.min()

    #################################################
    ################################################# Kolmogrov-Smirnov test
    ## using means.
    vals= (bedsel.IN + bedsel.OUT) / 2
    props= list(vals / chrom_end)
    vals= list(vals)

    kt= kstest(props, 'uniform')

    ### power of Kolmogorov-smirnov test
    ###

    stat_store, count_store= draw_comp_ks(n= bedsel.shape[0], rangepos= [chrom_start,chrom_end],
                  	comp= 'uniform',nrep= 10000,
                    geno_bornes= [chrom_start,chrom_end],Steps= ksteps,window_size=kwindow)

    cthist, ctidx= np.histogram(count_store,bins= len(set(count_store)), density= True)
    ct_cdf= [1 - (sum(cthist[:x]) / sum(cthist)) for x in range(1,len(cthist)+1)]

    #
    stat_store= pd.DataFrame(stat_store,columns= ['KS','PVAL','DF','DMEAN','DSTD'])
    stat_store.head()

    pvalktest= kstest(list(stat_store.PVAL),'uniform')

    power_here= stat_store[stat_store.PVAL < alpha].shape[0] / stat_store.shape[0]
    report= 'power (KS | H0): {} for alpha = {}'.format(stat_store[stat_store.PVAL < alpha].shape[0] / stat_store.shape[0], alpha)

    kt_out.append([chrom,bedsel.shape[0],*list(kt),power_here])


kt_out= np.array(kt_out)
kt_out= pd.DataFrame(kt_out, columns= ["CHROM","NINV","KS", "PVAL","FDR"])

kt_out.to_csv(out_dir + 'KStest.txt',index= False, sep= "\t")



