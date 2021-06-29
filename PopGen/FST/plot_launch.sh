#!/bin/bash
#SBATCH --time=5:00
#SBATCH --nodes=1
#SBATCH --partition=debug


module purge 
module load python/3.6.2

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ./env

python -u plot_fst.py
