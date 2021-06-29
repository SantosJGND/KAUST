#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --job-name=PLS
#SBATCH --mem=4GB

R/3.6.0/gnu-6.4.0

Rscript PCafindPLS.r
