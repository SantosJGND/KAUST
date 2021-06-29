#!/bin/bash
#SBATCH --time=5:00
#SBATCH --nodes=1
#SBATCH --partition=debug
SECONDS=0

module load vg/1.28.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /ibex/scratch/santj0a/GenomeGraphs/environments/env_odgi

graph=$1

echo $graph
homedir=`dirname $graph`

################################
######## VG

echo "view"
vg view -F -v $graph.gfa > $graph.vg

echo "mod"
vg mod -X 256 $graph.vg >$graph"_mod.vg"

echo "move"
mv $graph"_mod.vg" $graph.vg

echo "index"
vg index -x $graph.xg $graph.vg
echo "prune"
vg prune $graph.vg > $graph.pruned.vg

echo "index2"
vg index -g $graph.gcsa -Z 500 $graph.pruned.vg
rm -f $graph.pruned.vg

## SNARLS and PATHS
vg snarls -e -r $homedir"/traversals.pb" $graph.xg > $homedir"/snarls.pb"

vg find -I -x $graph.xg > $homedir"/paths.txt"

for path in `cat $homedir"/paths.txt"`; do
echo $path
vg deconstruct -p $path -r $homedir"/snarls.pb" -e -a $graph.xg > $homedir"/"$path.vcf
awk '(substr($0,1,1)!="#") {print $1,$2,$3,length($4),length($5)}' $homedir"/"$path.vcf > $homedir"/"$path"_SV.txt"
done

rm $homedir"/"*.vcf
################################
######### ODGI

odgi_dir=$homedir"/odgi/"
mkdir -p $odgi_dir

time odgi build -g $graph.gfa -o - -p | odgi sort -i - -o $graph.dg -b -w

# make a visualization
odgi viz -i $graph.dg -x 4000 -y 800 -L 0 -X 150 -P 20 -R -o $odgi_dir"graph.dg.png"
odgi layout -i $graph.dg -o $odgi_dir"graph2D.dg.svg"

echo "paths"
time odgi paths -i $graph.dg -L > $odgi_dir"paths_list.txt"

time odgi paths -i $graph.dg -H > $odgi_dir"haplotypes.txt"

echo "graph matrix"
time odgi matrix -i $graph.dg > $odgi_dir"graph.mat"


#####################################


ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

echo $ELAPSED

