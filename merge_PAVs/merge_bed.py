import sys
import itertools as it
#######################################################################
#######################################################################
from vector_process import *
from merge_tools import *
import pandas as pd
import numpy as np

import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--bedlist', type=str, default= 'PAV.bed')
parser.add_argument('--names', type=str, default= '')
parser.add_argument('--dist', type=int, default= 20)
parser.add_argument('--pavl', type=int, default= 50)
parser.add_argument('--mummer', type= str, default= '')

parser.add_argument('--id', type= str, default= 'Azucena')
parser.add_argument('-o', type= str, default= './')

args = parser.parse_args()

#################
flist= args.bedlist.split(',')

snames= args.names.split(',')
if len(snames) < len(flist):
        snames= ['s{}'.format(x) for x in range(len(snames))]

#################
pav_dict= {
        "INS": ["Ins","INS"],
        "DEL": ["DEL","Del"]
}

#################

mg_dict= import_split(flist,snames,id=args.id,pavs= pav_dict)


#bedins= pd.read_csv(bedinput, sep= "\t", header= None)
#bedins.columns= ["chrom", "IN", "OUT", "SV", "acc", "soft"]

query="soft"
dist_threshold= args.dist

for pavt, bedins in mg_dict.items():

        accessions= bedins[query].unique()
        accessions= list(accessions)
        #
        invbed= bedins.copy().reset_index(drop= True)
        invbed.columns= ["chrom","IN","OUT","SV","acc", "soft"]
        ### summary before merging 
        soft_summ= []
        for soft in invbed.soft.unique():
            soft_summ.append([soft,invbed[invbed.soft == soft].shape[0]])
        soft_summ= pd.DataFrame(soft_summ, columns= ["soft","total"])
        soft_summ.to_csv(args.o + pavt + "_{}_d{}.soft_count.bed".format('-'.join(snames),args.dist), index= False, sep= "\t")
        print("hello {}".format(invbed.shape))
        td= filter_info(invbed, accessions,query= query,dist_threshold= args.dist)

        td.to_csv(args.o + pavt + '_{}_d{}.merged.bed'.format('-'.join(snames),args.dist), index= False, sep= "\t")
        #
        mg_dict[pavt]= td

################################
################################
insl=mg_dict['INS']
dels= mg_dict['DEL']

print(insl.head())

PAV_azu=[]
for chrom in td.chrom.unique():

    pchrom= pd.concat((insl[insl.chrom == chrom],dels[dels.chrom == chrom])).sort_values('IN').reset_index(drop= True)
    PAV_azu.append(pchrom)

PAV_azu= pd.concat(tuple(PAV_azu)).sort_values('chrom').reset_index(drop= True)

print(PAV_azu.head())

#############
############# Subset Pandas Array
PAV_azu["L"]= pd.to_numeric(PAV_azu["L"])
pav_subset= PAV_azu[PAV_azu["L"] > 50].copy().reset_index(drop=True)

#############
pavarray= 1 - (pav_subset.iloc[:,4:].to_numpy() == "NA")
pavarray= np.array(pavarray,dtype= int)

print(pavarray.shape)
print(pavarray[:5])
################################
print(PAV_azu.shape)
softs= list(PAV_azu.columns[4:])
stats_comp(pavarray, fig=args.o + "{}_".format(args.dist), show= True,cols= softs)

percent= True

pav_stats= plot_bars(pavarray,snames,percent= percent, fig= args.o + "{}_".format(args.dist))
pav_stats.to_csv(args.o + 'summary_{}_d{}.tsv'.format('-'.join(snames),args.dist), index= False, sep= "\t")
################################
PAV_azu.to_csv(args.o + 'merged_{}_d{}.merged.bed'.format('-'.join(snames),args.dist), index= False, sep= "\t")
################################

if args.mummer:
	#### read mummer file
	mummer=pd.read_csv(args.mummer,sep="\t")
	mummer.columns= ["ref_start", "ref_end", "qry_start", "qry_end", "ref_len", "qry_len",
                  "identiy", "ref_tag","qry_tag"]

	mummer= mummer.sort_values('ref_start').reset_index(drop=True)
	
	pavdf= PAV_azu.iloc[:,:3].copy()
	pavdf.columns= ["chrom","start","end"]
	pavdf.chrom= ["Chr{}".format(x.zfill(2)) for x in pavdf.chrom]
	alt_coords= pav_adjust(pavdf, mummer, disth= 80, id=args.id, svh= args.o, pavl= args.pavl)
	print(alt_coords.head())
	print(alt_coords.shape)
	print(PAV_azu.shape)
	npc= pd.concat((PAV_azu,alt_coords),axis= 1)
	
	npc.to_csv(args.o + "mummer_assess" + '_{}_d{}.merged.bed'.format('-'.join(snames),args.dist), index= False, sep= "\t")


