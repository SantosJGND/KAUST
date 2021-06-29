#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --job-name=simeSMC
#SBATCH --mem=4GB


source /home/santj0a/miniconda3/etc/profile.d/conda.sh
conda activate /home/santj0a/Projects/SLiM_ABC/env_ms

module load R/3.6.0/gnu-6.4.0
module load tabix/0.2.6
module load bcftools/1.9/gnu-6.4.0

################################ DIR & PARAMS
vcf_sim=$1
project_name=`basename $vcf_sim`
project_name=${project_name%".vcf.gz"}

homedir=`pwd`

subdir=$project_name"/"
mkdir $subdir

inds_file=$homedir"/total_simSamples.txt"

Ninds=5
repeats=10

cd $subdir
################################### PREP Data
###################################
for rep in `seq 1 $repeats`; do

shuf -n $Ninds $inds_file > sample.txt

avail=""
##
for ind in `cat sample.txt`; do
bcftools view -s $ind ../$vcf_sim > sample_$ind.vcf
bgzip -f sample_$ind.vcf
tabix -p vcf sample_$ind.vcf.gz
avail=$avail" "sample_$ind.vcf.gz
done
##

python ../generate_multihetsep.py $avail > rep$rep"_multihetsep.txt"

done
rm sample_i*

#################################### RUN eSMC
####################################
Nhaps=$(($Ninds * 2))
tag=`basename $vcf_sim`
tag=${tag%".vcf"}

echo $tag

Rscript --vanilla ../run_eSMC.sh $repeats $Nhaps $tag

module purge
conda deactivate
