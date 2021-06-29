#!/bin/bash

#conda init
source /home/santj0a/miniconda3/etc/profile.d/conda.sh
echo `conda list | grep packages`
conda activate
conda activate /home/santj0a/Projects/SLiM_ABC/env_ms

echo `pwd`
#module load python/3.6.2

rescale=10

echo `conda list | grep packages`
python -u process_est.py $1 --temprec $2 --rescale $rescale

conda deactivate 

./slim -m -s $RANDOM instance_recipe.slim

