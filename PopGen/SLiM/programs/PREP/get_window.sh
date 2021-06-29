#!/bin/bash
#SBATCH --time=25:00
#SBATCH --nodes=1
#SBATCH --partition=debug

module purge
module load plink/5.3
module load angsd/1.0

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/santj0a/Projects/SLiM_ABC/env_dadi

##################### DIRECTORIES

project_dir=$1
scratchdir=$2
project_name=$3

subdir="VCF"

##################### FILES

plink_file="/home/santj0a/RICE_data/3K_plink/core_3K_nonCoding"

## remember:  samtools faidx fasta
fasta_ref="/home/santj0a/RICE_data/refseq/ncbi-genomes-2020-10-31/GCF_001433935.1_IRGSP-1.0_genomic.fna"
## gene ranges:
genranges="/home/santj0a/RICE_data/gff/GOpA_plink.txt"

temp_ids=$project_dir"/IDs/temp_IDs.txt"
subtrop_ids=$project_dir"/IDs/subtrop_IDs.txt"

###################### ARGUMENTS
samp_n=20
reps=25
maf=0.01
geno=0.2
window_size=100000

###################### GO
mkdir $subdir

##
################################### Create ind to group files by sampling
ind_file=$subdir"/ind_file.txt"
subprune=$subdir"/subprune"
ind_assign=$subdir"/ind_gp.txt"
tempfile=$subdir"/temp.txt"

vcf_file=$subdir"/gp"

#shuf -n $samp_n $temp_ids > $subdir"/gp1.txt"
#shuf -n $samp_n $subtrop_ids > $subdir"/gp2.txt"

cp $temp_ids $subdir"/gp1.txt"
cp $subtrop_ids $subdir"/gp2.txt"
Ntemp=`cat $temp_ids | wc -l`
Nsub=`cat $subtrop_ids | wc -l`

cat $subdir"/gp1.txt" > $ind_file
echo "" >> $ind_file
cat $subdir"/gp2.txt" >> $ind_file


#################################### Get random window
python window_collate.py --input $genranges --lkeep $window_size
## extract
plink --bfile $plink_file --keep $ind_file \
--maf $maf --geno $geno \
--extract range window_select.txt \
--recode vcf --out $vcf_file


################################### process ID names in files (underscore in 3K plink data);
sed 's/\t/_/g' $subdir"/gp1.txt" > $tempfile
python name_process.py $tempfile --tag gp1
cat $tempfile > $ind_assign

echo "" >> $ind_assign

sed 's/\t/_/g' $subdir"/gp2.txt" > $tempfile
python name_process.py $tempfile --tag gp2
cat $tempfile >> $ind_assign

################################### EASY SFS
###################################
python -u easySFS.py -i $vcf_file".vcf" -f -o $subdir"/out" \
-p $ind_assign --prefix $project_name --proj $Ntemp","$Nsub -a

cp $subdir"/out/fastsimcoal2/"$project_name"_jointMAFpop1_0.obs" ./
rm -r $subdir"/out"

################################### SUMMARY STATS
###################################
python fold_sfs.py $project_name"_jointMAFpop1_0.obs" OBStemp.txt \
--vcf $vcf_file".vcf" --inds $ind_assign

mv OBStemp.txt $project_name"_jointMAFpop1_0.obs"

###################################################################
#rm plink*

conda deactivate

