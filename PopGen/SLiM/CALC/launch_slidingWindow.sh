#!/bin/bash
#SBATCH --nodes=2
#SBATCH --time=3:00:00
#SBATCH --job-name=SIM_calc
#SBATCH --mem=8GB

source /home/santj0a/miniconda3/etc/profile.d/conda.sh
conda activate /ibex/scratch/santj0a/Projects/SLiM/libseq_env

ID=CrashSel
wind=100
stepsize=100

python statistics_slidingwindow_pylibseq_SingExon_osg.py -folder ../sims/ -simID $ID -winSize $wind -stepSize $stepsize

conda deactivate

