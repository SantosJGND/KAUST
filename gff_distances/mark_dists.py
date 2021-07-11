import pandas as pd 
import numpy as np


import argparse

parser = argparse.ArgumentParser(description='Process some integers.')

parser.add_argument('--gff1', type= str, default= "")
parser.add_argument('--gff2', type= str, default= "")
parser.add_argument('--f1', type= str, default= "gene")
parser.add_argument('--f2', type= str, default= "LTR")
parser.add_argument('-o', type= str, default= "summary.txt")

args = parser.parse_args()


genes= pd.read_csv(args.gff1, sep= "\t", header= None, dtype= str, comment= '#')
genes.columns= ["contig","m1","pred","start","end","q","s","t","idt"]
genes= genes[genes['pred'].str.contains(args.f1)].reset_index(drop= True)
genes.idt= [x.split(";")[0].split("=")[1] for x in genes.idt]

tes= pd.read_csv(args.gff2, sep= "\t", header= None, comment= '#', dtype= str)
tes.columns= ["contig","m1","pred","start","end","q","s","t","idt"]
tes= tes[tes['pred'].str.contains(args.f2)].reset_index(drop= True)
tes.idt= [x.split(";")[0].split("=")[1] for x in tes.idt]


summary_dump=[]

for chrom in genes.contig.unique():
	genc= genes[genes.contig == chrom].reset_index(drop= True)
	tec= tes[tes.contig == chrom].reset_index(drop= True)
	
	genpos= np.array(genc[["start","end"]], dtype= int)
	tepos= np.array(tec[["start","end"]], dtype= int)
		
	for ix in range(genc.shape[0]):
		start,end= genpos[ix,:]
		lpe= end - start
		
		if len(tepos) == 0:
			summ= [chrom, start, end, genc.idt[ix],"NA", "NA", 0, "NA"]
			summary_dump.append(summ)
			continue
		ovx= (tepos[:,0] < end) & (tepos[:,1] > start)
		overlap= tec[ovx]
		ovid= ",".join(overlap.idt)
		
		sdists= np.min(tepos - start, axis= 1)
		edists= np.min(tepos - end, axis= 1)
		dists= np.array([sdists, edists], dtype= float).T
		dists[dists < 0 ] = +np.Inf
		dists= np.min(dists,axis= 1)
		
		dists[ovx]= 0
		distmin= np.argmin(dists)
		summ= [chrom, start, end, genc.idt[ix],tec.idt[distmin], dists[distmin], len(overlap), ["NA",ovid][int(len(ovid) > 0)]]
		#print(summ)
		summary_dump.append(summ)
	
	
summary_dump= np.array(summary_dump, dtype= str)
summary_dump= pd.DataFrame(summary_dump)
summary_dump.columns= ["contig","start","end", "idt","closest","cdist","overlap","Overlap_id"]

summary_dump.to_csv(args.o,index= False, header= True, sep= "\t")


