#!/bin/bash

module purge
module load python/3.6.2


source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/x_garciaj/Projects/SLiM_ABC/tempsubtrop.v3/env


ind_assign=sim_ind_gp.txt
project_name=tempsubtrop.v3
samp_n=20

python -u easySFS.py -i $1 -f -o "out" \
        -p $ind_assign --prefix $project_name --proj $samp_n","$samp_n -a

cp "out/fastsimcoal2/"$project_name"_jointMAFpop1_0.obs" ./temp_obs.txt
rm -r out/

python fold_sfs.py ./temp_obs.txt
