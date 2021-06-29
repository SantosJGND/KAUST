#!/bin/bash

module load tabix/0.2.6

vcf=/ibex/scratch/projects/c2016/SANTOSJ_DIR/People/INVOICE/Alice/Ocoarctata_pop/vcf_input/Ocoarctata_b1_noRep_nomiss.vcf.gz
repeat_gff="/ibex/scratch/projects/c2016/SANTOSJ_DIR/People/INVOICE/Alice/Ocoarctata_pop/vcf_input/OcoaRS1.PZYT01.fa.out.gff"


outdir="/ibex/scratch/projects/c2016/SANTOSJ_DIR/Data/Ocoar_calls/"

zgrep "^#\|ChrK" $vcf | bgzip > $outdir"OcoarctataChrK_b1_noRep_nomiss.vcf.gz"
zgrep "^#\|ChrL" $vcf | bgzip > $outdir"OcoarctataChrL_b1_noRep_nomiss.vcf.gz"

echo -e "chrom\tstart\tend" > $outdir"repeat_mask.bed"
cut -f1,4,5 $repeat_gff | grep -v "^#" > $outdir"repeat_mask.bed"

module purge
