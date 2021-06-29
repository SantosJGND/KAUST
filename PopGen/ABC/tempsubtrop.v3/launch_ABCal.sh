#!/bin/bash
#SBATCH --time=5:00
#SBATCH --nodes=1
#SBATCH --partition=debug

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ./env

progdir_prep="/home/x_garciaj/Projects/SLiM_ABC/programs/PREP/"
progdir_sim=/home/x_garciaj/Projects/SLiM_ABC/programs/SIM/
progdir_summ=/home/x_garciaj/Projects/SLiM_ABC/programs/SUMM/

homedir=/home/x_garciaj/Projects/SLiM_ABC/tempsubtrop.v3.1

rec_dir=$homedir'/recipe/'
obs_dir=$homedir'/observed/'

##
input_file=tempsub_ABCcali.input

##
parfile=$rec_dir"tempsubtrop.v3.slim"
interface=$rec_dir"AB-Sinterface.txt"
estfile=$rec_dir"tempsubtrop.v3.est"
sim_ids=$rec_dir"sim_ind_gp.txt"
##
obsfile=$obs_dir"tempsubtrop.v3_jointMAFpop1_0.obs"

##
window_func=$progdir_prep"get_window.sh"
window_func2=$progdir_prep"window_collate.py"
window_func3=$progdir_prep"name_process.py"
##
simulationProg1=$progdir_sim"launch_SLiM.sh"
simulationProg2=$progdir_sim"slim"
int_prog=$progdir_sim"process_est.py"

##
SSProgram1=$progdir_summ"fold_sfs.sh"
SSProgram2=$progdir_summ"fold_sfs.py"
SSProgram3=$progdir_summ"easySFS.py"


###
##
scratchFolder="scratch"$1
mkdir $scratchFolder

homefolder=`pwd`

#copy files
cp $window_func $window_func2 $window_func3 $scratchFolder/
cp $input_file $interface $int_prog $parfile $estfile $obsfile $sim_ids $scratchFolder/
cp ABCsampler $simulationProg1 $simulationProg2 $SSProgram1 $SSProgram2 $SSProgram3 $scratchFolder/
#go on the node and launch ABCsampler1.0
cd $scratchFolder

chmod +x ABCsampler $simulationProgram $SSProgram get_window.sh
## launch
./get_window.sh $homefolder $scratchFolder

./ABCsampler $input_file addToSeed=1
#copy results back
cp *output*.txt *.log $homefolder/
