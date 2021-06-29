import numpy as np 
import os
import itertools as it
from scipy.stats import kstest,poisson 
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from sklearn.cluster import AgglomerativeClustering

def ins_labels(bedsel,inv_tab= ["IN","OUT"],dist_threshold= 1000.0):
    """

    """
    inv_arrays= [np.array(bedsel[y],dtype= int) for y in inv_tab]
    inv_labels= []

    for idx in range(len(inv_tab)):
        inv_distances= inv_arrays[idx]

        clustering = AgglomerativeClustering(n_clusters= None,distance_threshold= dist_threshold)

        clusters= clustering.fit(inv_distances.reshape(-1,1))
        inv_labels.append(clusters.labels_)


    inv_labels= np.array(inv_labels).T
    
    return inv_labels

##
def merge_labs(label_array):
    '''
    merge array of label columns.
    '''
    
    unique_now= label_array[0,:]
    nlabs= []
    labi= 0
    
    for idx in range(label_array.shape[0]):
        comp_now= label_array[idx,:]
        comparison= comp_now == unique_now
        comparison= sum(comparison)
        
        if comparison == 0:
            labi += 1
            nlabs.append(labi)
        else:
            nlabs.append(labi)
        
        unique_now= comp_now
    
    return nlabs

def bed_collapse(bed_array,clcol= 'INV_merge',start_cols= ['chrom','IN'],
                 end_cols= ['OUT'],func_do= [np.mean, np.mean]):
    '''
    collapse array based on labels column.
    '''
    #print(bed_array.head)
    inv_clusts= list(bed_array[clcol])
    inv_dict= {z: [x for x in range(bed_array.shape[0]) if inv_clusts[x] == z] for z in set(inv_clusts)}
    clust_set= set(inv_clusts)
    new_array= []
    new_cols= []
    new_cols.extend(start_cols)
    new_cols.extend(end_cols)
    
    new_cols.extend(['L', "idx", "trace"])
    #new_cols.append('idx')
    for cl,g in inv_dict.items():
        new_row= []
        minmax= []
        for col in start_cols:
            ends= [bed_array.iloc[x][col] for x in g]
            
            if col == 'IN':
                ends= int(func_do[0](ends))
                minmax.append(ends)
            else:
                ends= ends[0]
            
            new_row.append(ends)
        
        for col in end_cols:
            ends= [bed_array.iloc[x][col] for x in g]
            
            if col == 'OUT':
                ends= int(func_do[1](ends))
                minmax.append(ends)
            else:
                ends= ends[0]
            
            new_row.append(ends)
        
        if len(minmax)== 2:
            le= minmax[1] - minmax[0]
            
            new_row.append(int(le))
        
        clumped=[]
        for clp in g:
        	nc=[bed_array.iloc[clp]["chrom"],bed_array.iloc[clp]["IN"], bed_array.iloc[clp]["OUT"]]
        	nc= [str(x) for x in nc]
        	nc= "-".join(nc)
        	clumped.append(nc)

        new_row.append(".".join([str(x) for x in g]))
        new_row.append(';'.join([str(x) for x in clumped]))
        new_array.append(new_row)
    
    new_array= np.array(new_array)
    
    new_array= pd.DataFrame(new_array,columns= new_cols)
    return new_array


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

parser.add_argument('--kwindow', type=int,
                    default=2e5)

parser.add_argument('--ksteps', type=int,
                    default=1e5)

parser.add_argument('--alpha', type=float,
                    default=1e-2)

parser.add_argument('--cmet', type=str,
                    default="fdr_bh")

parser.add_argument('--tdist', type=float,
                    default=15000.0)

parser.add_argument('--centro', type="str",
                    default="Nipp.Centro.txt")

args = parser.parse_args()



kwindow= args.kwindow
ksteps= args.ksteps
alpha= args.alpha
corr_method= args.cmet
dist_threshold= args.tdist


file_dir= 'C:/Users/floyd/Desktop/KAUST/Projects/INV_distribution/INV_windows_27-03-2021/'
file= file_dir + 'filtered_HQ_INV.bed'

###
if args.centro:
    centromeres= pd.read_csv(centrom_file,header= None,sep= '\t')
    centromeres.columns= ['ID','chrom','IN','OUT','L']

###
invbed= pd.read_csv(file, header=None, sep= '\t')
invbed.columns= ['chrom','IN','OUT','L']

## replace names in bed column
#invbed.chrom= invbed.chrom.str.split('.',n= 1, expand= True)[1]

invbed.head()


out_dir= '/'.join(file.split('/')[:-1]) + "/"+ file.split('/')[-1].split('.')[0] + '_{}_{}_{}_{}'.format(kwindow,alpha,corr_method,args.tdist) + "/"


#######################
#######################
pdout= []
kt_out= []
bed_filtered= []
for chrom in invbed.chrom.unique():

    bedsel= invbed.loc[invbed.chrom==chrom].reset_index(drop= True)
    bedsel= bedsel.sort_values(by="IN")
    proxy_bed= bedsel.copy()

    INV_labels= ins_labels(proxy_bed, dist_threshold= dist_threshold)

    proxy_bed["INV_merge"]= merge_labs(INV_labels)

    bedsel= bed_collapse(proxy_bed)
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

    ### plot KS stats
    figname= out_dir + 'KS_stats/{}_KSsummary.pdf'.format(chrom)
    os.makedirs(os.path.dirname(figname), exist_ok=True)

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))

    axes[0,0].hist(stat_store.DMEAN, label= 'mean')
    axes[0,0].legend()
    axes[0,1].hist(stat_store.DSTD,label= 'std')
    axes[0,1].legend()

    axes[1,0].scatter(stat_store.DMEAN*1e7,stat_store.DSTD*1e7)
    plt.setp(axes[1,0], xlabel='mean',ylabel= 'std')

    axes[1,1].hist(stat_store.PVAL)
    plt.setp(axes[1,1], xlabel='pval')

    plt.savefig(figname)
    plt.close()
    ####################### ################################
    ###################### #################################
    ###################### #################################

    db_russex= bedsel.copy()
    db_russex.columns= ['chrom','POS','OUT','L','idx','trace']
    #
    db_russex.head()

    kn= db_russex.POS.shape[0]

    ## DETERMINE RATE
    kn= db_russex.POS.unique().shape[0]
    lamb= kn * kwindow /(chrom_end - chrom_start)

    # WINDOWS
    Windows,Out= geno_Lwind_split(db_russex, geno_bornes= [chrom_start,chrom_end],Steps= ksteps,window_size=kwindow)

    # OBS CDF
    winds= sorted(Windows.keys())
    lss= [len(Windows[x]) for x in winds]
    lss= np.array(lss)

    pt= [cthist[0] / sum(cthist)] + ct_cdf + [0] * (max(lss) - len(cthist))
    ########### window K cdf
    ###########
    fig_cdf_name= out_dir + 'Wcdf/{}_Wcdf.pdf'.format(chrom)
    os.makedirs(os.path.dirname(fig_cdf_name), exist_ok=True)

    xs= list(range(0,len(pt)))

    ##
    pdftest = 1 - poisson.cdf(xs,lamb)

    proxypd= np.array(xs,dtype=float)

    pdftest[proxypd <= lamb]= poisson.cdf(proxypd[proxypd <= lamb],lamb)
    pdftest[proxypd > lamb]= 1 - poisson.cdf(proxypd[proxypd > lamb],lamb)
    #

    plt.figure(figsize= (5,5))
    plt.plot(xs,pdftest)

    plt.title("nINV cdf; wL: {} rate: {}".format(kwindow,lamb))
    plt.ylabel("P")
    plt.xlabel("K")
    plt.ylim(0,1.1)

    plt.savefig(fig_cdf_name)
    plt.close()

    ############
    # MULTIPLE TEST CORRECTION

    from statsmodels.stats import multitest
    pdftest= np.array(pt)[lss]
    pdf_bool, pdf_correct, alpha_S, alpha_B= multitest.multipletests(pvals=pdftest, alpha= alpha, method=corr_method)

    alpha_sel= alpha_S

    #########
    ######### INDEX PVALS 
    add_invs= True
    width, height= [18,5]
    #
    pdf_dict= np.array(pdf_bool,dtype= int)

    pdf_dict= {z: [x for x in range(len(winds)) if pdf_dict[x] == z] for z in [0,1]}

    ########## PLOT
    ##########
    figname= out_dir+ 'dist/{}_dist.pdf'.format(chrom)

    os.makedirs(os.path.dirname(figname), exist_ok=True)

    fig, ax = plt.subplots(1,figsize=(width, height))

    ## CENTROMERE
    if args.centro:
        ##
        centro= centromeres[centromeres.chrom == chrom]
        centro_l= centro.OUT - centro.IN
        centro_l= list(centro_l)[0]

        centri= list(centro.IN)[0]
        errorboxes= []
        rect = Rectangle((centri, min(pdftest)), centro_l, max(pdftest)-min(pdftest))
        errorboxes.append(rect)

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes)

    # Add collection to axes
    ax.add_collection(pc)

    ## INVERSIONS
    if add_invs:
        for INV in bedsel.IN:
            plt.vlines(int(INV),min(pdftest),max(pdftest)-min(pdftest),alpha=0.3)

    ## WINDOW CORRECTED PVAL
    col_th= ['black','red']
    for z,g in pdf_dict.items():
        plt.scatter(np.array(winds)[g], pdftest[g],color= col_th[z])

    plt.title('{} window analysis'.format(chrom))
    plt.xlabel('positions')
    plt.ylabel('window pmf')
    plt.savefig(figname)
    plt.close()

    if 1 in pdf_dict.keys():
        INV_select= list(it.chain(*[Windows[winds[x]] for x in pdf_dict[1]]))
        INV_select= list(set(INV_select))

        INV_select= db_russex.iloc[INV_select,:]
        original_idx= list(INV_select.idx)
        
        original_idx= it.chain(*[[int(x) for x in y.split('.')] for y in original_idx])
        original_idx= sorted(original_idx)
        
        out_select= proxy_bed.iloc[original_idx,:]
        pdout.append(out_select)


if len(pdout):
	pdout= pd.concat(pdout)
	pdout.to_csv(out_dir + 'selected.txt',index= False, sep= "\t")

bed_filtered= pd.concat(bed_filtered)
bed_filtered.to_csv(out_dir + "INV_concat.bed", index= False, sep= "\t")


kt_out= np.array(kt_out)
kt_out= pd.DataFrame(kt_out, columns= ["CHROM","NINV","KS", "PVAL","FDR"])

kt_out.to_csv(out_dir + 'KStest.txt',index= False, sep= "\t")



