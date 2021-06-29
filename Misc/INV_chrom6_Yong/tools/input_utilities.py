
import re
import numpy as np
import pandas as pd 
import gzip

######################################### INPUT


def read_focus(index_file):
    indxs = []
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        indxs.append(line[0])
    
    Input.close()
    
    return indxs


def VCF_read_filter(sim, sim_dir= './',chrom= '1',haps_extract= False, scale_genSize= False,
    collapsed= True,min_size= 5, samp= [30,30,5], stepup= 'increment', outemp= './', ploidy= 2, return_sub= False,
    indfile= 'ind_assignments.txt', fasta_template= 'chr{}_{}.fa.gz', diffs= False, bases= 'ACGT', ksize= 3):
    '''
    read vcf file. filter snps. 
    '''
    vcf_dir= sim_dir + sim + '/'
    vcf_file= vcf_dir + sim + '_' + 'chr' + chrom + '.vcf.gz'
    
    print(sim)

    #### read vcf
    genotype, summary, Names= read_vcf_allel(vcf_file,haps_extract= haps_extract)   
    
    print(genotype.shape)
    print(summary.shape)

    if len(genotype) == 0:
        return {}, {}, {}, {}, {}, {}
        
    ## read fasta
    fasta_file= vcf_dir + fasta_template.format(chrom,sim)

    with gzip.open(fasta_file,'r') as f:
        lines= f.readlines()
        lines= [x.decode() for x in lines]

    refseq= lines[1].strip()

    #### scale mutation counts. 
    L= len(refseq)
    scale= 1
    if scale_genSize:
        if genome_size > 1:
            scale= int(genome_size/L)

    ### subset window (for application to full data sets) if necessary.
    ### Not used here, wstart and wend set to min and max positions respectively.
    positions= [int(x) for x in summary.POS]
    wstart= int(min(positions))-1
    wend= int(max(positions))
    
    Wlen= wend - wstart
    
    genotype_parse= [x for x in range(summary.shape[0]) if int(summary.POS[x])-1 >= wstart and int(summary.POS[x])-1 <= wend]
    Window= genotype[:,genotype_parse]
    subset_summary= summary.loc[genotype_parse,:].reset_index()
    
    ### get mutation-type by SNP matrix, 
    ### filter SNPs if necessary (see vcf_muts_matrix_v1)
    t0= time.time()
    mut_matrix, flag_reverse, flag_remove= vcf_muts_matrix_v1(refseq,subset_summary,start= wstart,end= wend,ksize= ksize,
        bases=bases, collapsed= collapsed)
    
    print('mut_matrix shape: {}'.format(mut_matrix.shape))
    retain= [x for x in range(Window.shape[1]) if x not in flag_remove]
    Window= Window[:,retain]
    subset_summary= subset_summary.loc[retain,:].reset_index()

    t1= time.time()
    time_mut= t1 - t0

    if diffs:
    	sim_start= sim.split('.')[-1]
    	diff_snps= read_diffs(sim,diff_dir= vcf_dir, start= int(sim_start))

    	summary_diff= [x for x in range(subset_summary.shape[0]) if subset_summary.POS[x] in diff_snps.keys()]

    	flag_reverse.extend(summary_diff)
    	flag_reverse= list(set(flag_reverse))
    
    
    if flag_reverse:
        Window[:,flag_reverse]= ploidy - Window[:,flag_reverse]
    #
    if return_sub:
        return Window, mut_matrix, scale, refseq, subset_summary
    else:
        return Window, mut_matrix, scale




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

