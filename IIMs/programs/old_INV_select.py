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
submat.to_csv(home + "INV_sampSel_kept{}.txt".format(min_inds), sep= "\t")

keep_invs= np.random.randint(0,submat.shape[0],size= keep)
keep_invs= sorted(keep_invs)
submat= submat.iloc[keep_invs,:].reset_index()

inversions= submat.iloc[:,5:]

####
####


files_content= []

idx= 0
row= inversions.iloc[idx]
print(genlist)

###
control= inversions.columns[row.isnull()]
control= list(control)
files_control= genlist.code.isin(list(control))
files_control= genlist.irgsp[files_control]
files_control= list(files_control)
#files_control= [seqs_home + x for x in files_control]

###
outside= inversions.columns[~row.isnull()]
outside_files= genlist.code.isin(outside)
outside_files= genlist.irgsp[outside_files]
outside_files= list(outside_files)
#outside_files= [seqs_home + x for x in outside_files]

####
chrom= submat.CHROM[0]
chrom_here= chrom_IRGSP.default==chrom
chrom_here= chrom_IRGSP.irgsp[chrom_here]
chrom_irgsp= list(chrom_here)[0]

start= submat.Start[0] - bornes
start= int(start)
end= submat.End[0] + bornes
end= int(end)

length_inv= submat.Length[0] + bornes * 2
length_inv= int(length_inv)

#### PRINT

outdir="base_info/"

with open(outdir + "ref_fasta.txt", 'w') as fp:
	line= '\n'.join(["\t".join([x,x]) for x in outside_files])
	fp.write(line)

with open(outdir + "control_fasta.txt", 'w') as fp:
	line= '\n'.join(["\t".join([x,x]) for x in files_control])
	fp.write(line)

with open(outdir + "region_length.txt", 'w') as fp:
	fp.write(str(length_inv))

with open(outdir + "region.txt", 'w') as fp:
	line= '\t'.join([chrom,str(start),str(end),"inv"])
	fp.write(line)
