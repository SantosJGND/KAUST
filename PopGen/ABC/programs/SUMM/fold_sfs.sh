#!/bin/bash

module purge
module load python/3.6.2


source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/santj0a/Projects/SLiM_ABC/env_dadi

project_name=tempsubtrop.v3

################################################
ind_obs=VCF/ind_gp.txt
ind_assign=sim_ind_gp.txt
>$ind_assign

samp1=`cat VCF/gp1.txt | wc -l`
samp2=`cat VCF/gp2.txt | wc -l`

d=0
for i in `seq 0 $(($samp1-1))`; do
echo -e i$d"\tgp1" >> $ind_assign
d=$(($d+1))
done

for i in `seq 0 $(($samp2-1))`; do
echo -e i$d"\tgp2" >> $ind_assign
d=$(($d+1))
done


##################################################

python -u easySFS.py -i $1 -f -o "out" \
-p $ind_assign --prefix $project_name --proj $samp1","$samp2 -a

cp "out/fastsimcoal2/"$project_name"_jointMAFpop1_0.obs" ./temp_obs.txt
rm -r out/

python fold_sfs.py ./temp_obs.txt --vcf $1 --inds $2

conda deactivate
