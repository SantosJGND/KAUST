#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=3:00:00
#SBATCH --job-name=mergeS
#SBATCH --mem=45GB

module load bedtools
source ~/miniconda3/etc/profile.d/conda.sh

conda activate /ibex/scratch/projects/c2016/SANTOSJ_DIR/conda/environments/python_env

dir_graph=/ibex/scratch/projects/c2016/SANTOSJ_DIR/People/Yong/graphM16_PAVformat/
dir_trace=/ibex/scratch/projects/c2016/SANTOSJ_DIR/People/Yong/graphM16_PAVformat/
dir_svim=/ibex/scratch/santj0a/GenomeGraphs/VCF_work/VCF_wg/PAV/
dir_asm=/ibex/scratch/projects/c2016/SANTOSJ_DIR/People/Yong/svim_asm/assemblies/bed_format/
dir_snif=/ibex/scratch/projects/c2016/SANTOSJ_DIR/Projects/PAVs/sniffles/
dir_vulcan=/ibex/scratch/projects/c2060/SANTOSJ_DIR/Projects/SV_calling/Vulcan/

d=2
flank=15

outdir=single_merge_asmmg/
mkdir -p $outdir

outsumm=$outdir"summary.txt"

tractfile=Azucena_tract.bed

for tractfile in `cat names_order.txt`; do

name=${tractfile%"_tract.bed"}
subname=${name#Os}

### which files correspond to this accession (different names)
tractfile=$dir_graph$tractfile
asmf=$dir_asm`ls $dir_asm | grep $subname -`
svimf=$dir_svim"genome1.genome"$d".svim_pavSumm.txt"
sniff=$dir_snif"genome1.genome"$d".NGMLR.sniffles.bed"
vulcan=$dir_vulcan"genome1.genome5.vulcan.sniffles.bed"

coord_file=/ibex/scratch/projects/c2016/Yong/S06_Joao_NucmerResults/genome1.genome"$d".filter.coord

echo $tractfile $asmf

echo $name
echo $outdir 

python -u merge_bed.py --bedlist $tractfile","$asmf --names "mg,asm" \
--dist $flank --pavl 50 -o $outdir$name"_" --id $name \
#--mummer $coord_file

d=$((d+1))

done
