#!/bin/bash
#SBATCH --nodes=2
#SBATCH --time=150:00:00
#SBATCH --job-name=AB-Sr
#SBATCH --mem=8GB

./ABCsampler $1 addToSeed=1

#copy results back
cp *output*.txt *.log $homefolder/

