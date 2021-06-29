import pandas as pd
import numpy as np
import itertools as it
import collections
import re
import scipy

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import estimate_bandwidth
from sklearn.cluster import MeanShift, estimate_bandwidth

from structure_tools.Modules_tools import return_fsts2

from structure_tools.AMOVA_func import AMOVA_FM42, amova_cofactor


def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


def geno_subset_random(genotype, summary, RG_info, ID_col,subset_col,Names,include= [],code= {},others= 'admx',Sn= 500, Sm= 10000):

    ### Subset to acceptable range of accessions x markers.

    Present= [x for x in range(len(Names)) if Names[x] in list(RG_info[ID_col])]
    if include:
        presinc= [x for x in include if x in Names]

        presidx= [x for x in range(len(Names)) if Names[x] in presinc] 
        print('{}/{} of requested names present.'.format(len(presinc),len(include)))

    if len(Present) < genotype.shape[0]:
        '{} IDs missing'.format(genotype.shape[0] - len(Present))

    Nsample= sorted(np.random.choice(Present,Sn,replace= False))
    
    if include:
        Nsample.extend(presidx)
        Nsample= sorted(list(set(Nsample)))

    Msample= sorted(np.random.choice(list(range(genotype.shape[1])),Sm,replace= False))

    ###
    gen_sample= np.array(genotype[Nsample,:])

    gen_sample= np.array(gen_sample[:,Msample])

    subsummary= summary.loc[Msample,:]

    subsummary= subsummary.reset_index()

    Names_select= [Names[x] for x in Nsample]

    ###

    print('gen_sample shape: {}, {}'.format(len(Nsample),len(Msample)))

    ###

    Name_idx= [list(RG_info[ID_col]).index(x) for x in Names_select]
    code_vec= [RG_info[subset_col][x] for x in Name_idx]

    #code_vec= [code_vector[x] for x in Nsample]
    code_vec= [[x,others][int(x not in code.keys())] for x in code_vec]

    code_vec= [code[x] for x in code_vec]

    code_lib= {
        z:[x for x in range(len(code_vec)) if code_vec[x] == z] for z in list(set(code_vec))
    }

    return gen_sample, subsummary, code_vec, code_lib, Nsample, Msample



def read_geno_nanum(filename, row_info= 6,header_info= 9,phased= False):

    info_summ= {}
    info_save= list(range(row_info))

    header_len= header_info
    summary= []

    Miss= recursively_default_dict()
    
    Input= open(filename)

    genotype= []
    d= 0

    for line in Input:    
        line= line.strip()

        if d in info_save:

            line= ''.join(filter(lambda ch: ch not in "#", line))
            line= line.split('=')
            info_summ[line[0]] = ''.join(line[1:])
            d += 1
            continue

        if d == (len(info_save)):
            print(info_summ)
            line= ''.join(filter(lambda ch: ch not in "#", line))
            line= line.split()

            columns= line[:header_len]
            Names= line[header_len:]

            d += 1
            continue

        if d > (len(info_save)):
            line= line.split()
            seq= []

            info= line[:header_len]
            chrom= re.search(r'\d+', line[0]).group()
            info[0]= chrom

            summary.append(info)

            for ind in range(header_len,len(line)):
                locus= line[ind]
                #print(locus)
                alleles= locus.split(':')[0]
                #print(alleles)
                if '.' in alleles:
                    alleles= ''.join([[x,'0'][int(x == '.')] for x in list(alleles)])
                alleles= list(map(int, re.findall(r'\d+', alleles)))
                if len(alleles) != 2:
                    print(alleles)
                if phased:
                    seq.extend(alleles)
                else:
                    seq.append(sum(alleles))

            genotype.append(seq)
            d += 1

    Input.close()

    summary= np.array(summary)
    summary= pd.DataFrame(summary,columns= columns)
    genotype= np.array(genotype).T

    return genotype, summary, Names


def simple_read_vcf(filename,row_info= 5,header_info= 9,phased= False):

    Input= open(filename)

    info_summ= {}
    info_save= list(range(row_info))

    phased= False
    header_len= header_info
    summary= []

    Miss= recursively_default_dict()

    genotype= []
    d= 0

    for line in Input:    
        line= line.strip()

        if d in info_save:

            line= ''.join(filter(lambda ch: ch not in "#", line))
            line= line.split('=')
            info_summ[line[0]] = ''.join(line[1:])
            d += 1
            continue
        
        if d == (len(info_save)):
            line= ''.join(filter(lambda ch: ch not in "#", line))
            line= line.split()

            columns= line[:header_len]

            Fam= {
                line[x]: x for x in range(header_len,len(line))
            }

            d += 1
            continue

        if d > (len(info_save)):
            line= line.split()
            seq= []

            info= line[:header_len]
            chrom= re.search(r'\d+', line[0]).group()
            info[0]= chrom

            summary.append(info)
            
            for ind in range(header_len,len(line)):
                locus= line[ind]
                #print(locus)
                alleles= locus.split(':')[0]
                #print(alleles)
                if '.' in alleles:
                    alleles= ''.join([[x,'0'][int(x == '.')] for x in list(alleles)])
                alleles= list(map(int, re.findall(r'\d+', alleles)))
                if len(alleles) != 2:
                    alleles
                if phased:
                    seq.extend(alleles)
                else:
                    seq.append(sum(alleles))

            genotype.append(seq)
            d += 1

    Input.close()

    summary= np.array(summary)
    summary= pd.DataFrame(summary,columns= columns)
    genotype= np.array(genotype).T

    return genotype, summary, info_save


def check_densities(vector_lib_2,N):
    
    who= np.random.choice(list(range(vector_lib_2.shape[0])),N,replace= False)
    
    freqs= []
    
    for pop in who:
        
        freq_vector= vector_lib_2[pop,:]

        X_plot = np.linspace(0, 1, 100)

        kde = KernelDensity(kernel='gaussian', bandwidth=0.01).fit(np.array(freq_vector).reshape(-1,1))

        log_dens= kde.score_samples(X_plot.reshape(-1,1))
                        
        freqs.append(np.exp(log_dens))
    
    freqs= np.array(freqs)
    
    return freqs


def geno_window_split(genotype,summary,Steps= 25,window_size=100):

    window_starts= list(np.arange(0,genotype.shape[1],Steps))
    if window_starts[-1] != genotype.shape[1]:
        window_starts.append(genotype.shape[1])

    Windows= recursively_default_dict()
    Out= recursively_default_dict()

    lengths_winds= []

    for splyt in range(len(window_starts) - 1):
        IN= window_starts[splyt]
        OUT= window_starts[splyt] + window_size

        if OUT > genotype.shape[1]:
            OUT= genotype.shape[1] - 1
        range_window= [IN,OUT]

        lengths_winds.append(OUT-IN)

        chrom= int(summary.CHROM[range_window[0]])

        start= int(summary.POS[range_window[0]])
        end= int(summary.POS[range_window[1]])

        Windows[chrom][start]= genotype[:,range_window[0]:range_window[1]]
        Out[chrom][start]= end

        if OUT-IN < window_size:
            break

    return Windows, Out


def window_fst_sup(Windows,ref_labels,labels1,Chr= 1,ncomp= 4,range_sample= [],rand_sample= 0):
    
    kde_class_labels= labels1
    kde_label_dict= {
        z:[x for x in range(len(kde_class_labels)) if kde_class_labels[x] == z] for z in list(set(kde_class_labels))
    }
    
    if rand_sample:
        sample= rand_sample
        sample_range= [0,sample]
        Freq_extract= {
            Chr:{
                bl:Windows[Chr][bl] for bl in np.random.choice(list(Windows[Chr].keys()),sample,replace= True)
            }
        }

    if range_sample:
        sample_range= range_sample
        Freq_extract= {
            Chr:{
                bl:Windows[Chr][bl] for bl in list(sorted(Windows[Chr].keys()))[sample_range[0]:sample_range[1]]
            }
        }
    
    sim_fst= []
    
    for c in Freq_extract[Chr].keys():
        Sequences= Windows[Chr][c]

        if Sequences.shape[1] <= 3:
            Results[Chr][c] = [0,0]
            print('hi')
            continue

        Sequences= np.nan_to_num(Sequences)

        pca = PCA(n_components=ncomp, whiten=False,svd_solver='randomized').fit(Sequences)
        data = pca.transform(Sequences)
        Ngps= len(ref_labels)
        these_freqs= []
        
        for hill in ref_labels:
            cl_seqs= Sequences[kde_label_dict[hill],:]

            freq_vector= [float(x) / (cl_seqs.shape[0] * 2) for x in np.sum(cl_seqs,axis= 0)]
            these_freqs.append(freq_vector)
            
        Pairwise= return_fsts2(np.array(these_freqs))
        sim_fst.append(list(Pairwise.fst))
    
    ###
        
    return sim_fst


