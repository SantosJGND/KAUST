#!/bin/bash

###
### MODULES
module load minimap2/2.17
module load samtools/1.8

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /ibex/scratch/projects/c2016/SANTOSJ_DIR/conda/environments/skbio_env
###
### FILES

bam_dir=/ibex/scratch/projects/c2016/Yong/S01_Joao_10samplesBamFrom3K/
refseq=/ibex/scratch/projects/c2016/SANTOSJ_DIR/Data/refseq/GCF_001433935.1_IRGSP-1.0_genomic.fna.gz

invbed=inv_intervals.bed

###
### PARAMS
qualt=40


###
### RUN

#for inv in `cat $invbed`; do 

inv="Chr01:4723836-4725649"

python refsubset.py --refseq $refseq --reg $inv -o inv_fasta.fa

#for bam in $bam_dir; 
#if [[ ! $acc *bai ]]; then

bam=$bam_dir"XI-3A_ERS470559.sorted.bam"

acc=${bam%.*}
acc=`basename $acc`

samtools view -h -q $qualt $bam $inv | samtools bam2fq - > $acc.fq

minimap2 -xsr inv_fasta.fa $name.fq > $acc.paf


module purge
