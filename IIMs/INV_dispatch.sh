#!/bin/bash
#SBATCH --time=5:00
#SBATCH --nodes=1
#SBATCH --partition=debug

module purge
module load plink/5.3
module load vcftools/0.1.17

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /ibex/scratch/santj0a/GenomeGraphs/environments/python_env
refplink=/home/santj0a/RICE_data/3K_plink/NB_final_snp

run_name="run"$2

homedir=`pwd`
program_dir=$homedir"/programs/"

mkdir -p $run_name"/"
cd $run_name

cp $program_dir"parse_snps.py" ./
cp $program_dir"INV_select.py" ./

mkdir -p base_info

#################
#################
keep=$1

inv_matrix=$homedir"/prelim_data/INV_sampSel_kept3.txt"

python -u INV_select.py --matrix $inv_matrix --keep $keep

cp base_info/ref_fasta.txt base_info/combined_samp.txt
cat base_info/control_fasta.txt >> base_info/combined_samp.txt


############################################### SUBSET AND GROUP FREQUENCY
###############################################

plink --bfile $refplink --keep base_info/combined_samp.txt \
--extract range base_info/region.txt \
--geno 0.2 \
--maf 0.01 \
--make-bed --out inv_subset_one

plink --bfile inv_subset_one \
--geno 0.2 \
--maf 0.01 \
--make-bed --out inv_subset

plink --bfile inv_subset \
--recode vcf \
--out inv_subset

plink --bfile inv_subset --keep base_info/ref_fasta.txt \
--freq

sed 's/ \+ /\t/g' plink.frq | sed 's#^\t##g' | sed 's#^ ##g' | cut -f1-5  > ref.freq
sed -i -e 's#MAF#ALT1#g' ref.freq
head ref.freq
plink --bfile inv_subset --keep base_info/control_fasta.txt \
--freq

mv plink.frq ctrl.freq
sed -i -e 's#MAF#ALT2#g' ctrl.freq

##################################
################################## VCFTOOLS
grep "#" inv_subset.vcf > temp.vcf
grep -v "#" inv_subset.vcf | sed 's/^/Chr/g' >> temp.vcf
mv temp.vcf inv_subset.vcf

sed -e 's#\t#_#g' base_info/control_fasta.txt > base_info/control_list.txt
sed -e 's#\t#_#g' base_info/ref_fasta.txt > base_info/ref_list.txt

vcftools --vcf inv_subset.vcf  --weir-fst-pop base_info/control_list.txt --weir-fst-pop base_info/ref_list.txt --out control_vs_ref

sed 's/ \+ /\t/g' ctrl.freq | sed 's#^\t##g' | sed 's#^ ##g' | cut -f5 | paste ref.freq - > temp1

cut -f3 control_vs_ref.weir.fst | paste temp1 - > summary.txt
sed -i -e '1 ! s/^/Chr/g' summary.txt

#################
#################

python -u parse_snps.py --summ summary.txt -o snplist

if [ -s snplist ]; then 

plink --bfile $refplink \
--extract snplist \
--geno 0.3 \
--maf 0.001 \
--make-bed --out inv_3K
#--extract range base_info/region.txt \

plink --bfile inv_3K --r2
mv plink.ld selected.ld
rm plink.*
fi

#################
rm inv* parse_snps.py ctrl.freq ref.freq temp1 control*

module purge
