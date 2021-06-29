import numpy as np
import pandas as pd

seqs_home="/home/santj0a/RICE_data/refseq/"
home="/ibex/scratch/santj0a/GenomeGraphs/VG/INV_subdir/IMMs/prelim_data/"
inv_matrix= home + "INV_matrix.txt"
genomes_list= home + "genomes_list.txt"
irgsp_chroms= home + "IRGSP_chroms.txt"

chrom_IRGSP=pd.read_csv(irgsp_chroms,sep= "\t")
genlist= pd.read_csv(genomes_list,sep="\t")
invmat= pd.read_csv(inv_matrix, sep="\t")

###
###
bornes= 5e3
keep= 1
min_inds= 3

available= genlist.copy()[~genlist["file"].isnull()]
keep_cols= list(invmat.columns[:4]) + list(np.array(sorted(available.code),dtype= str))[1:]

submat= invmat.copy()[keep_cols]

inversions= submat.iloc[:,5:]

inv_keep = (inversions.isnull().sum(axis=1) < (inversions.shape[1] - min_inds)) & (inversions.isnull().sum(axis=1) > min_inds)

submat= submat[inv_keep]
submat.to_csv(home + "INV_sampSel_kept{}.txt".format(min_inds), sep= "\t", index= False)

