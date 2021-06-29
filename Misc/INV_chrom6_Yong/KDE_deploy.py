
from tools.pretty_tools import pretty_input
from tools.input_utilities import read_focus
import itertools as it

import argparse

parser = argparse.ArgumentParser()


parser.add_argument('--vcf', type=str,
                    default='.')

parser.add_argument('--info', type=str,
                    default='3K_info.txt')

parser.add_argument('--IDcol', type=str,
                    default='ID')

parser.add_argument('--refs', type= str,
					default= '0,1,2', help= "comma del ref group labels in info table / IDcol column.")  

parser.add_argument('--chr', type=int,
                    default=6)

parser.add_argument('--ws', type=int,
                    default=150, help= 'window size')

parser.add_argument('--step', type=int,
                    default=15, help= 'step size')

parser.add_argument('--vc', type=int,
                    default=5, help= 'number of pca comps to keep locally.')

parser.add_argument('--rep', type=int,
                    default=20, help= 'number of repeat numbers ')

parser.add_argument('--refcol', type=str,
                    default='Initial_subpop')

parser.add_argument("--focus",type= str,
                    help = "reference accessions indexes in genofile.")

parser.add_argument('--ploidy', type=int,
                    default=2)

parser.add_argument('--fixed', action= 'store_true',
                    default= False)

parser.add_argument('--haps', action= 'store_true',
                    default= False)

parser.add_argument('--same', action= 'store_false',
                    default= True)

parser.add_argument('--sn', type=int,
                    default=0, help= "number of haps to sample per group, 0 = use ref dict.")

parser.add_argument('--out_dir', type=str,
                    default='./')

parser.add_argument("--MSprint",action= "store_true",help = "if given prints cluster stats.")


parser.add_argument("--id",type= str,default= '1',help = "Give your analysis an ID. default is set to integer 2")


args = parser.parse_args()



#####

from dict_code import (
    code
) 

ID_col= args.IDcol
subset_col= args.refcol
CHR= args.chr
window_size= args.ws
Steps= args.step
fixed= args.fixed
ploidy= args.ploidy

ref_labels= [int(x) for x in args.refs.split(',')]

#####
import os
start= args.id

Home= args.out_dir

print('writting analysis id:{0} to directory {1}'.format(args.id,Home))

filename= Home + "Blocks_Request_st"+str(start)+"_CHR" + str(CHR).zfill(2) + ".txt"
print(filename)
os.makedirs(os.path.dirname(filename), exist_ok=True)
###########
###########

genotype, summary, Names, RG_info= pretty_input(args.vcf,info_file= args.info,ID_col= ID_col,
                                                haps_extract= args.haps, ploidy= ploidy)   


#####
#####
print(Names[:10])
if args.focus:
    Names_select = read_focus(args.focus)
else:
    Names_select = list(Names)


if args.haps and args.focus:
    Names_select= [[y + '_{}'.format(x) for x in range(ploidy)] for y in Names_select]
    Names_select= list(it.chain(*Names_select))


Names_idx= [Names.index(x) for x in Names_select]

#####
#####
from tools.pretty_pack import window_parse

Windows, Out= window_parse(genotype, summary, fixed= fixed, window_size= window_size, Steps= Steps)

print('number of chromosomes: {}'.format(len(Windows)))
print('number of windows: {}'.format(sum([len(Windows[x].keys()) for x in Windows.keys()])))

## get profiles
from tools.KDE_tools import KDE_window_profiles

Topo, var_comp_store, Construct= KDE_window_profiles(genotype,Windows,Names_idx,
                    RG_info, ID_col,subset_col,Names,
                    ref_labels, Chr= CHR,n_comps= args.vc, exclude= [], gp_focus= [],
                    repn= args.rep, code= code,others= 'admx',Sn= args.sn, same= args.same,
                    ms_comp= args.MSprint)

##

###############################
###############################
import os
start= args.id

Home= args.out_dir

print('writting analysis id:{0} to directory {1}'.format(args.id,Home))

filename= Home + "Blocks_Request_st"+str(start)+"_CHR" + str(CHR).zfill(2) + ".txt"
print(filename)
os.makedirs(os.path.dirname(filename), exist_ok=True)
Output = open(filename,"w")

Output.write("CHR\tIn\tOut\tRef\t")


for var in Names_select:
    Output.write(var + "\t")

Output.write("\n")

for block in range(len(Topo[CHR][0])):
    for ref in Topo[CHR].keys():
        Output.write(str(CHR) + "\t")
        Output.write(str(int(Points[block])) + "\t")
        Output.write(str(int(Out[CHR][Points[block]])) + "\t")
        Output.write(str(ref) + '\t')
        for ass in range(len(Topo[CHR][ref][block])):
            Output.write(str(Topo[CHR][ref][block][ass]) + "\t")
        Output.write("\n")

Output.close()

#
#
#

if args.MSprint == True:

    filename= Home + 'Blocks_profiles_st'+str(start)+'_CHR'+ str(CHR).zfill(2)+ '.txt'
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    Output= open(filename,'w')

    Output.write('CHR\tIN\tcluster\t')

    for var in Names_select:
        Output.write(var + "\t")

    Output.write("\n")

    for prf in Construct[CHR].keys():
        for cl in Construct[CHR][prf].keys():
            Output.write(str(CHR) + "\t")
            Output.write(str(int(prf)) + '\t')
            Output.write(str(cl) + '\t')
            Output.write('\t'.join([str(round(x,5)) for x in Construct[CHR][prf][cl]]))
            Output.write('\n')

    Output.close()

