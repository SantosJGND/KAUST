#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ./env


gff=test_data/hapsolo.gff3
tempd=`basename $gff`".temp"

min_size=10

awk '$5 - $4 > '$min_size' {print $1 "\t" $4 "\t" $5 "\t" $5-$4 }' $gff > $tempd

python -u Distribution_simple.py --input $tempd

rm $tempd
