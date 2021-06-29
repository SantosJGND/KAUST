
import pandas as pd
import numpy as np
import itertools as it

import matplotlib.pyplot as plt


import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--input', type=str,
                    default="GOpA_plink.txt")

parser.add_argument('--lkeep', type=int,
                    default=2e5)

parser.add_argument('--out', type=str,
                    default="window_select.txt")

args = parser.parse_args()


ranges_db= args.input

lkeep= args.lkeep


### input
ranges_db= pd.read_csv(ranges_db,sep='\t',header= None)
ranges_db.columns= ["CHROM","IN","OUT","ID"]

chrom= np.random.choice(ranges_db.CHROM.unique(),size= 1)[0]

bedsel= ranges_db.loc[ranges_db.CHROM == chrom].copy().reset_index()
bedsel= bedsel.sort_values(by=['IN'])


### Merge overlapping genes
###
last= bedsel.loc[0]["OUT"]

merged_db= [[bedsel.loc[0]['CHROM'],0,bedsel.loc[0]['IN'],bedsel.loc[0]["IN"], bedsel.loc[0]["ID"]]]

for idx in range(1,bedsel.shape[0]):
    if bedsel.loc[idx]["IN"] > last:
        
        merged_db.append([bedsel.loc[idx]["CHROM"], last, bedsel.loc[idx]["IN"], bedsel.loc[idx]["IN"] - last,bedsel.loc[idx]["ID"]])
        last= bedsel.loc[idx]["OUT"]

merged_db= np.array(merged_db)
merged_db= pd.DataFrame(merged_db,columns= ["CHROM","IN","OUT","L","ID"])
merged_db["prop"]= np.array(merged_db.L,dtype= int) / sum(np.array(merged_db.L,dtype= int))
merged_db.head()


### merge windows between genes
###


start= np.random.choice(list(range(merged_db.shape[0])), size= 1, p= np.array(merged_db.prop, dtype= float))

si= start[0]
keep= [list(merged_db.loc[si])]

cl= merged_db.loc[start[0]]["L"]
cl= int(cl)

direct= [-1,1]
idx= [-1,1]

borders= [0, merged_db.shape[0]]
toggle= 0

while cl < lkeep:
    #
    dt= idx[toggle]
    
    nidx= si + dt
    
    if nidx <= borders[0]:
        idx[0]= 1
        direct[0]= 1
        toggle= 1 - toggle
        continue
        
    if nidx >= borders[1]:
        idx[1]= -1
        direct[1]= -1
        toggle= 1 - toggle
        continue
    
    nrow= merged_db.loc[nidx].copy()
    
    nl= int(nrow.loc["L"])
    
    nsum= nl + cl
    
    if nsum > lkeep:
        ###
        ###
        sL= nsum - lkeep
        nl= int(nl - sL)
        
        nrow.loc["OUT"]= int(nrow.loc["IN"]) + nl
        nrow.loc["L"]= nl

    
    keep.append(list(nrow))
    cl += nl
    
    idx[toggle] += direct[toggle]
    
    toggle= 1 - toggle


    
keep= np.array(keep)
keep= pd.DataFrame(keep, columns= ["CHROM","IN","OUT","L","ID","P"])
keep= keep.sort_values(by=['IN'])


keep.loc[:,["CHROM","IN","OUT","ID"]].to_csv(args.out,sep= "\t",header= None,index= False)


