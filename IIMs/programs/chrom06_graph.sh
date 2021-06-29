#!/bin/bash
#SBATCH --nodes=2
#SBATCH --time=6:00:00
#SBATCH --job-name=BeagleDemTest
#SBATCH --mem=9GB

module load minimap2/2.17
module load seqwish/0.6
module load edyeet/0.1

seq_use=/ibex/scratch/santj0a/GenomeGraphs/VG/chrom6_edyeet/chr06.fa.gz

#pggb_home=/home/santj0a/GenomeGraphs/VG/pggb/
#cp $pggb_home"pggb" ./

seq_use=$1
similarity=70
overlap=2000

##
echo $seq_use
echo $similarity
echo $overlap
##

dir_seq=`dirname $seq_use`"/"

./pggb_minimap -i $seq_use -s $overlap -p $similarity -a $similarity -n 10 -t 16 -v -l -S -K 18
