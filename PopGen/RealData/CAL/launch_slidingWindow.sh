#!/bin/bash
#SBATCH --nodes=2
#SBATCH --time=3:00:00
#SBATCH --job-name=SIM_calc
#SBATCH --mem=8GB

source /home/santj0a/miniconda3/etc/profile.d/conda.sh
conda activate /ibex/scratch/santj0a/Projects/SLiM/libseq_env

ID=CrashSel
wind=1000
stepsize=1000

python statistics_slidingwindow_pylibseq_SingExon_osg.py -folder ../temp_windows/ -simID $ID -winSize $wind -stepSize $stepsize

conda deactivate

