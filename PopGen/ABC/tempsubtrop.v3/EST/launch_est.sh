#!/bin/bash
#SBATCH --nodes=2
#SBATCH --time=6:00:00
#SBATCH --job-name=ABCest
#SBATCH --mem=7GB

./ABCestimator estimator.input
