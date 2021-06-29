#!/bin/bash
#SBATCH --nodes=1
#SBATCH -c 25
#SBATCH --time=15:50:00
#SBATCH --job-name=JellyfishGsE
#SBATCH --mem=150GB

module purge
module load R/3.6.0/gnu-6.4.0

##################################### K-mer size

kmer=17

##################################### PREP
## Input and output dirs
dir=/ibex/scratch/projects/c2016/SANTOSJ_DIR/GSest_demo/paeruginosa-reads/filters/
input=$dir*".filter.fastq"

odir="GenomeSizeEst/"$kmer"_clean/"
ofile=$odir$kmer"mer_out"

##################################### Genomescope

## Genome scope dir
gsdir=/ibex/scratch/projects/c2016/SANTOSJ_DIR/TOOLS/genomescope2.0/

Rscript $gsdir"genomescope.R" \
-i $ofile.histo \
-k $kmer \
-o $odir \
-m 1000000 \
-p 1

