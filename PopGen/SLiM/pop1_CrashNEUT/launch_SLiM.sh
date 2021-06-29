#!/bin/bash
#SBATCH --nodes=2
#SBATCH --time=25:00:00
#SBATCH --job-name=CrashNEU
#SBATCH --mem=8GB

#conda init
source /home/santj0a/miniconda3/etc/profile.d/conda.sh
echo `conda list | grep packages`
conda activate
conda activate /home/santj0a/Projects/SLiM_ABC/env_ms

echo `pwd`
#module load python/3.6.2

rescale=10

####################
tag=CrashBGS
repn=200

####################
sim_dir="sims/"

new_dir=`basename $1`
new_dir=${new_dir%"_interface.txt"}
echo $new_dir
simdir=sims/

for rep in `seq 1 $repn`; do 

tag=$new_dir"_"$rep

echo `conda list | grep packages`
python -u process_est.py $1 --temprec $2 --rescale $rescale --out $simdir --tag $tag

time ./slim -m -s $RANDOM $simdir$tag.slim

rm ANC_trees/$tag.trees

done

conda deactivate
