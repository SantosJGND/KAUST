
import numpy as np

from tools.parse_tools import (window_split_idx,geno_Lwind_split)



def window_parse(genotype, summary, fixed= True, Steps= 10,window_size= 120,
                    min_snp= 11,
                geno_bornes= []):
    ##### 
    if fixed:

        Windows, Out= window_split_idx(genotype,
                                        summary,
                                        Steps= Steps,
                                        window_size=window_size)

    else:

        Windows, Out= geno_Lwind_split(summary,geno_bornes= [],
                                     Steps= Steps, window_size= window_size)

        Windows, Out= wind_compress(Windows,Out,min_snp= min_snp)

        Windows= {Chr: Windows}
        Out= {Chr: Out}

    return Windows, Out


from tools.Sim_ideogram_tools import (
	compress_ideo, class_and_ideo, return_ideogram
	)

from tools.KDE_tools import KDE_window_profiles

def quick_kdeClass(genotype, summary,Names_idx,RG_info, Names, Chr= 6,subset= [], ref_labels= [0,1,2], n_comps= 4, ID_col= "IRIS_ID",
                  fixed= True,window_size= 120, Steps= 10, Names_select= [],subset_col= 'Initial_subpop',code= {},
                  Outlier_threshold= 1e-4,Comparison_threshold= 4,Sn= 90, repn= 20,alt_col= [],
                  id_tag= 'kde', color_lookup= {}, square= False, start= 0, end= 0, clean= False):
    
    #####
    ## split into windows

    Windows, Out= window_parse(genotype, summary, fixed= fixed, window_size= window_size, Steps= Steps)

    print('number of chromosomes: {}'.format(len(Windows)))
    print('number of windows: {}'.format(sum([len(Windows[x].keys()) for x in Windows.keys()])))

    ## get profiles
    Windows_profiles, var_comp_store= KDE_window_profiles(genotype,Windows,Names_idx,
                        RG_info, ID_col,subset_col,Names,
                        ref_labels, Chr= Chr,n_comps= n_comps, exclude= Names_select, gp_focus= [],
                        repn= repn, code= code,others= 'admx',Sn= Sn, same= True, clean= clean)

    ##

    groups_plot= [0]

    Blocks, ideo_kde, chromosome_list = class_and_ideo(Windows_profiles,Out,[0]*len(Names_select),
                             Names_select= Names_select,
                             Comparison_threshold= Comparison_threshold,
                             Outlier_threshold= Outlier_threshold,
                             groups_plot=groups_plot,
                             colors= 'standard',
                             alt_col= alt_col)

    ###


    ID= 'kde{}_gp{}_w{}_N{}_M{}_th{}'.format(id_tag,'-'.join([str(x) for x in groups_plot]),Sn,len(Names_select),window_size,Comparison_threshold)  


    Fig_ideo= return_ideogram(ideo_kde,chromosome_list,ID,color_lookup= color_lookup,
                              height= 20,width= 40,yfont= 25,square= True, start= subset[0],end= subset[1],
                                 xticks= 2.5e5,xfont= 20)

