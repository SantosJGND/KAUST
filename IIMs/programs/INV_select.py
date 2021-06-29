import numpy as np
import pandas as pd

seqs_home="/home/santj0a/RICE_data/refseq/"
home="/ibex/scratch/santj0a/GenomeGraphs/VG/INV_subdir/IIMs/prelim_data/"
genomes_list=home + "genomes_list.txt"

irgsp_chroms= home + "IRGSP_chroms.txt"
chrom_IRGSP=pd.read_csv(irgsp_chroms,sep= "\t")
genlist= pd.read_csv(genomes_list,sep="\t")

###
###
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--matrix', default= 'INV_matrix.txt', type=str,
                    help='matrix of inversions')

parser.add_argument('--keep', type= int, default=0,
                    help='index of matrix to keep')

args = parser.parse_args()


matrix_read= args.matrix
bornes= 0
keep= args.keep

submat= pd.read_csv(matrix_read, sep= "\t")

submat= submat.iloc[[keep],:].reset_index()

inversions= submat.iloc[:,5:]

####
####
files_content= []

idx= 0
row= inversions.iloc[idx]
###
control= inversions.columns[row.isnull()]
control= list(control)
files_control= genlist.code.isin(list(control))
files_control= genlist.irgsp[files_control]
files_control= list(files_control) + ["CX140"]
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

with open(outdir + "control_fasta.txt", 'w') as fp:
	line= '\n'.join(["\t".join([x,x]) for x in outside_files])
	fp.write(line)

with open(outdir + "ref_fasta.txt", 'w') as fp:
	line= '\n'.join(["\t".join([x,x]) for x in files_control])
	fp.write(line)

with open(outdir + "region_length.txt", 'w') as fp:
	fp.write(str(length_inv))

with open(outdir + "region.txt", 'w') as fp:
	line= '\t'.join([chrom,str(start),str(end),"inv"])
	fp.write(line)
