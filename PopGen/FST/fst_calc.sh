#!/bin/bash
#SBATCH --time=5:00
#SBATCH --nodes=1
#SBATCH --partition=debug

module purge
module load vcftools/0.1.17
module load plink/5.3
module load python/3.6.2

project_dir=`pwd`"/"

maf=0.01

## source data
plink_file="/home/x_garciaj/RICE_data/3K_plink/core_3K_nonCoding"
ind_file=$project_dir"IDs/temp_subtrop_IDs.txt"

temp_ids=$project_dir"IDs/temp_IDs.txt"
subtrop_ids=$project_dir"IDs/subtrop_IDs.txt"

## temp data
focus_vcf="tempsub_dat"
gp1_sub=$project_dir"temp_ids.txt"
gp2_sub=$project_dir"subtrop_ids.txt"

## out data
focus_vcf="tempsub_dat"


###
### Process
# get single column IDs to pass to vcf-tools
cut -f1 $temp_ids > "temp.txt"
paste -d "_" temp.txt temp.txt > $gp1_sub
cut -f1 $subtrop_ids > "temp.txt"
paste -d "_" temp.txt temp.txt > $gp2_sub

rm temp.txt

###
### 
plink --bfile $plink_file \
--keep $ind_file \
--maf $maf \
--indep-pairwise 50 5 0.5 \
--recode vcf --out $focus_vcf

vcftools --vcf $focus_vcf".vcf" \
--weir-fst-pop $gp1_sub \
--weir-fst-pop $gp2_sub
