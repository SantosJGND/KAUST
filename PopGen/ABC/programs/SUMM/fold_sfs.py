import os
from os import path
import numpy as np
import pandas as pd
import itertools as it

from sklearn.metrics import pairwise_distances

import argparse

def get_fst(genotype_array, ploidy= 2):

	props= [np.sum(x,axis= 0) / (x.shape[0] * ploidy) for x in genotype_array]

	pairwise_fst= return_fsts2(np.array(props))

	return pairwise_fst.loc[0]["fst"]


def return_fsts2(freq_array):
    pops= range(freq_array.shape[0])
    H= {pop: [1-(freq_array[pop,x]**2 + (1 - freq_array[pop,x])**2) for x in range(freq_array.shape[1])] for pop in range(freq_array.shape[0])}
    Store= []

    for comb in it.combinations(H.keys(),2):
        P= [sum([freq_array[x,i] for x in comb]) / len(comb) for i in range(freq_array.shape[1])]
        HT= [2 * P[x] * (1 - P[x]) for x in range(len(P))]
        per_locus_fst= [[(HT[x] - np.mean([H[p][x] for p in comb])) / [HT[x],1][int(HT[x] == 0)],0][int(HT[x] == 0)] for x in range(len(P))]
        per_locus_fst= np.nan_to_num(per_locus_fst)
        Fst= np.mean(per_locus_fst)

        Store.append([comb,Fst])
    
    return pd.DataFrame(Store,columns= ['pops','fst'])



def PI_Taj(dataM):
    
    n= dataM.shape[0]
    pairDiff= pairwise_distances(dataM, metric='manhattan')
    mu= np.tril_indices(n)
    pair_Diff= sum(pairDiff[mu])
    
    pi= (2 / n / (n-1)) * pair_Diff
    
    return pi

## GGVE p. 62 eq. 2.40
def Watt_est(dataM, L= 0):
    
    Sn= np.sum(dataM,axis= 0)
    Sn= [x for x in Sn if x > 0]
    Sn= len(Sn)
    
    An= [1 / j for j in range(1,dataM.shape[0])]
    An= sum(An)
    
    Tw= Sn / An
    
    return Tw

## GGVE p. 62 eq. 2.42
def TajD(dataM):
    
    n= dataM.shape[0]
    
    Pit= PI_Taj(dataM)
    Watt= Watt_est(dataM)
    
    ##
    Sn= np.sum(dataM,axis= 0)
    Sn= [x for x in Sn if x >= 0]
    Sn= len(Sn)
    
    An= [1 / j for j in range(1,n)]
    An= sum(An)
    
    Bn= [1 / j**2 for j in range(1,n)]
    Bn= sum(Bn)
    
    E1= (n + 1) / (3 * An * (n-1)) - 1 / An**2
    E2= 2 * (n**2 + n + 3) / 9*n / (n-1)
    E2= E2 - (n+2)/ n / An + Bn / An**2
    
    E2= E2 / (An**2 + Bn)
    
    ##
    denom= E1 * Sn + E2 * Sn * (Sn - 1)
    denom= np.sqrt(denom)
    
    D= Pit - Watt
    D= D / denom
    
    return D


def popPA(dataM, pop_dict):
	totalpa= np.sum(dataM,axis= 0)
	SEGdict= {z: np.sum(dataM[g,:],axis= 0) for z,g in pop_dict.items()}
	PAdict= {z: [x for x in range(dataM.shape[1]) if g[x] == totalpa[x]] for z,g in SEGdict.items()}
	PAdict= {z: len(g) for z,g in PAdict.items()}
	
	return PAdict


os.environ["NUMEXPR_MAX_THREADS"]="272"
import allel

def read_vcf_allel(file_vcf,haps_extract= False,calldata= 'calldata/GT'):
    '''
    Use scikit allel to read vcf file. Organise variant information into summary pandas df. 
    '''
    geno1= []

    vcf_ori= allel.read_vcf(file_vcf)

    if not vcf_ori:
        print('file:')
        print(file_vcf)
        print('is empty.')

        return {}, {}, {}

    ### get genotype array
    geno= vcf_ori[calldata]
    
    mult_alt= []
    indel= []
    single= []

    ## Filter SNPs. append to single list what to 
    for idx in range(geno.shape[0]):
        ## eliminate +1 segregating mutations.
        if vcf_ori['variants/ALT'][idx][1]:
            gen_t= geno[idx]
            gen_t[gen_t > 1] = 0
            geno[idx]= gen_t
            ## or just jump them
            indel.append(idx)

        elif len(vcf_ori['variants/REF'][idx]) != 1 or len(vcf_ori['variants/ALT'][idx][0]) != 1:
            indel.append(idx)
        else:
            single.append(idx)

    if haps_extract:
        geno1= geno[:,:,0].T
        geno= geno[:,:,1].T
        geno= np.concatenate((geno,geno1),axis= 0)
    else:
        geno= allel.GenotypeArray(geno)
        geno= geno.to_n_alt().T

    ## setup summary

    column_names= ['CHROM','POS','ID','REF','ALT','QUAL','FILTER']

    alts= [vcf_ori['variants/ALT'][x][0] for x in range(geno.shape[1])]
    PASS= [['.','PASS'][int(vcf_ori['variants/FILTER_PASS'][x])] for x in range(geno.shape[1])]

    summary= [
        vcf_ori['variants/CHROM'],
        vcf_ori['variants/POS'],
        vcf_ori['variants/ID'],
        vcf_ori['variants/REF'],
        alts,
        vcf_ori['variants/QUAL'],
        PASS,

    ]

    summary= np.array(summary).T

    if len(indel):
        #
        geno= geno[:,single]
        if len(geno1):
            geno1= geno1[:,single]
        summary= summary[single,:]

    summary= pd.DataFrame(summary,columns= column_names)
    
    return geno, summary, vcf_ori['samples']


##################################################
##################################################



parser = argparse.ArgumentParser()

parser.add_argument("input",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read. Any number can be given.")

parser.add_argument("--vcf",type= str,default= '',
	help = "sequence length")

parser.add_argument("--inds",type= str,default= '',
	help = "inds to group")

parser.add_argument("--sfs", action= "store_true",default= False,
        help = "include sfs")


args = parser.parse_args()



vcf_file= args.input[0]

if len(args.input) == 1:
	out_file= "summary_stats_temp.txt"
else:
	out_file= args.input[1]


with open(vcf_file,'r') as fp:
	lines= fp.readlines()
	lines= [x.strip().split() for x in lines]

lines= lines[2:]
lines= np.array(lines,dtype= str)
lines= lines.T[1:]


if args.vcf:

	vcf_array= args.vcf

	genotype, summary, Names= read_vcf_allel(vcf_array,haps_extract= False)

	if args.inds:
		with open(args.inds,'r') as fp:
			inds_l= fp.readlines()
			inds_l= [x.strip().split() for x in inds_l]
			gps= list(set([x[1] for x in inds_l]))
			gp_dict= {z: [x for x in range(len(inds_l)) if inds_l[x][1] == z] for z in gps}

	else:
		gp_dict= {'gp0': list(range(genotype.shape[0]))}

	print("####")
	print(genotype.shape)
	datas= [genotype[gp_dict[z]] for z in sorted(gp_dict.keys())]
	
	##
	fst= get_fst(datas)
	##
	PA_dict= popPA(genotype,gp_dict)
	PAnames= sorted(PA_dict.keys())
	##
	stats_get= []

	for dat in datas:
		homm_stat= dat == 0
		homm_stat= np.sum(homm_stat) / np.prod(homm_stat.shape)
		stats_get.extend([homm_stat,PI_Taj(dat), TajD(dat)])


with open(out_file,'w') as fp:

	if args.vcf:

		for gp in sorted(gp_dict.keys()):
			head= [x.format(gp) for x in ['{}_H','{}_Tpi','{}_TajD']]
			head= '\t'.join(head)
			fp.write(head + '\t')

		fp.write('FST' + '\t')
		fp.write('\t'.join(['PA_{}'.format(x) for x in PAnames]))

	
	if args.sfs:
		for i in range(lines.shape[0]):
			for j in range(lines.shape[1]):
				fp.write('f{}_{}'.format(i,j) + '\t')

	fp.write('\n')

	if args.vcf:
		fp.write('\t'.join([str(x) for x in stats_get]) + '\t')
		fp.write(str(fst) + '\t')
		fp.write('\t'.join([str(PA_dict[x]) for x in PAnames]))

	if args.sfs:
		for i in range(lines.shape[0]):
			for j in range(lines.shape[1]):
				fp.write(str(lines[i,j]) + '\t')

print('SFS folded')


