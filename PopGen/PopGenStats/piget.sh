#SBATCH --nodes=1
#SBATCH --time=00:40:00
#SBATCH --job-name=vcfPI
#SBATCH --mem=8GB

module purge
module load vcftools/0.1.17
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /ibex/scratch/projects/c2016/SANTOSJ_DIR/conda/environments/python_env

vcf=/ibex/scratch/projects/c2016/SANTOSJ_DIR/Data/Ocoar_calls/OcoarctataChrL_b1_noRep_nomiss.vcf.gz
repeat_mask=/ibex/scratch/projects/c2016/SANTOSJ_DIR/Data/Ocoar_calls/repeat_mask.bed

name=`basename $vcf`
name=${name%.vcf.gz}

if [ ! -f NEstats.txt ]; then
echo -e "file\tPIm\tNEte\tNEadh" > NEstats.txt
fi

for keep in IDs/*; do
gpfile=`basename $keep`
gpname=${gpfile%.*}

vcftools --gzvcf $vcf --keep $keep --TajimaD 100000 --exclude-bed $repeat_mask --out $name"_"$gpname
vcftools --gzvcf $vcf --keep $keep --window-pi 100000 --exclude-bed $repeat_mask --out $name"_"$gpname

python draw_pi.py --pi $name"_"$gpname.windowed.pi
python draw_Taj.py --taj $name"_"$gpname.Tajima.D

done

rm *log
rm *.pi
rm *.D

conda deactivate
module purge

