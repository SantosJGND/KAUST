#!/bin/bash

vcf_dir="/mnt/d/GitHub/KAUST/Misc/INV_chrom6_Yong/vcf_data/"
vcf_file=$vcf_dir"chr6_phased_II.vcf.gz"


rg_file="/mnt/d/GitHub/KAUST/Misc/INV_chrom6_Yong/3K_info.txt"

echo $vcf_file

python -u KDE_deploy.py \
--vcf $vcf_file \
--info $rg_file \
--IDcol IRIS_ID \
--refs 0,1,2 \
--chr 6 \
--ws 150  \
--step 75 \
--rep 1 \
--MSprint \
--fixed \
--haps 
