#!/bin/bash
#SBATCH --time=5:00
#SBATCH --nodes=1
#SBATCH --partition=debug


homedir=/home/santj0a/Projects/SLiM_ABC/tempsubtrop.v2.1
scratchdir=/ibex/scratch/santj0a/Projects/SLiM/tempsubtrop.v2.1/

progdir_prep="/home/santj0a/Projects/SLiM_ABC/programs/PREP/"
progdir_sim=/home/santj0a/Projects/SLiM_ABC/programs/SIM/
progdir_summ=/home/santj0a/Projects/SLiM_ABC/programs/SUMM/

rec_dir=$homedir'/recipe/'
obs_dir=$homedir'/observed/'

##################################################################### Input
#####################################################################
input_file=tempsub_ABCcali.input

##
parfile=$rec_dir"tempsubtrop.v2.1.slim"
interface=$rec_dir"AB-SinterfaceR.txt"
estfile=$rec_dir"tempsubtrop.v2.1.est"
#sim_ids=$rec_dir"sim_ind_gp.txt"

##
window_func=$progdir_prep"get_window.sh"
window_func2=$progdir_prep"window_collate.py"
window_func3=$progdir_prep"name_process.py"
##
simulationProg1=$progdir_sim"launch_SLiM.sh"
simulationProg2=$progdir_sim"slim"
simulationProg3=$progdir_sim"launch_sampler.sh"
int_prog=$progdir_sim"process_est.py"

##
SSProgram1=$progdir_summ"fold_sfs.sh"
SSProgram2=$progdir_summ"fold_sfs.py"
SSProgram3=$progdir_summ"easySFS.py"


####################################################################### ENV setup
##
project_name=`basename $homedir`
scratchFolder=$scratchdir"scratch"$1
mkdir $scratchFolder

homefolder=`pwd`

#copy files
cp $window_func $window_func2 $window_func3 $scratchFolder/
cp $input_file $interface $int_prog $parfile $estfile $scratchFolder/
cp ABCsampler $simulationProg1 $simulationProg2 $simulationProg3 $SSProgram1 $SSProgram2 $SSProgram3 $scratchFolder/
#go on the node and launch ABCsampler1.0
cd $scratchFolder

chmod +x ABCsampler $simulationProgram $SSProgram get_window.sh

#########################################################################
######################################################################### launch
./get_window.sh $homedir $scratchFolder $project_name

sbatch launch_sampler.sh $input_file
