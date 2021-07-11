#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda activate /ibex/scratch/projects/c2016/SANTOSJ_DIR/conda/environments/python_env


## input files
home=/ibex/scratch/projects/c2016/For_João/

edta=O_ridleyi_canu_ordered_checked_GAP_fixed.fasta.mod.EDTA.intact.gff
edta=$home$edta
genes=round2.all.maker.noseq.gff
genes=$home$genes

## output
output=summary_dists.txt


###
python -u mark_dists.py --gff1 $genes --gff2 $edta -o $output
 


