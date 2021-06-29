#!/bin/bash
#SBATCH --nodes=1
#SBATCH -c 25
#SBATCH --time=15:50:00
#SBATCH --job-name=JellyfishGsE
#SBATCH --mem=150GB

module purge
module load jellyfish/2.2.10


##################################### K-mer size

kmer=17

##################################### PREP
## Input and output dirs
dir=/ibex/scratch/projects/c2016/SANTOSJ_DIR/GSest_demo/paeruginosa-reads/filters/
input=$dir*"trimmed.fastq"

odir="GenomeSizeEst/"$kmer"_clean/"
ofile=$odir$kmer"mer_out"

##
###################################### RUN
mkdir -p $odir

jellyfish count \
-t 25 \
-C \
-m $kmer \
-s 1000000000 \
-o $ofile \
--min-qual-char=? $input

jellyfish histo -o $ofile.histo $ofile --high=100000
