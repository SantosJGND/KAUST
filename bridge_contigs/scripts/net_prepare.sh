#!/bin/bash
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --time=0:15:00
#SBATCH --job-name=GnetHap
#SBATCH --mem=5GB

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /ibex/scratch/projects/c2016/SANTOSJ_DIR/conda/environments/python_env

pb_list=pb.fofn
pri_asm=../../quast/Cgileadensis.hifiasm.p_ctg.2e5_new.fa_1000_0.5988_0.6430to1.8910_0.2111_primary.fasta
#pri_asm=../../haplomerge/HaploMerger2_20180603/project_CG/genome_A_ref.fa
### STEP 1.

qual=$1
rlen=$2
minD=$3
rprop=$4

Op=$5
Ot=$6
Om=$7

minmap=150

gff=../../EDTA/Cgileadensis.hifiasm.p_ctg.2e5_new.fa_1000_0.5988_0.6430to1.8910_0.2111_primary.fasta.mod.EDTA.intact.gff3

edge_files=""

params=""
for i in $qual $rlen $minD $rprop $Op $Ot $Om; do params=$params$i"\t"; done


dir="parse_q"$qual"_l"$rlen"_m"$minD"_r"$rprop"_Op"$Op"_Ot"$Ot"_Om"$Om"/"
mkdir -p $dir


for i in `cat $pb_list`;

do
	#dir="parse_q"$qual"_l"$rlen"_m"$minD"_r"$rprop"/"
	#mkdir -p $dir
	name=`basename $i` ## 
	
	links=$name.map_merge.tsv
	name=$dir$name
	##
	cat $pri_asm | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > $dir"ref_contiglens.tsv"
	##
	awk '$7 > '$qual' && $4 > '$rlen' && $6 > '$minmap' && $6 < ($4 * '$rprop')' $links > $name.filtered.tsv
	
	awk '$5 < ($4 * '$rprop')' $name.filtered.tsv > $name.edgein
	awk '$5 > ($3 - ($4 * '$rprop'))' $name.filtered.tsv > $name.edgeout

	echo "edgins: "`cat $name.edgein | wc -l`
	echo "edgeouts: "`cat $name.edgeout | wc -l`
	
	join -o 1.1,1.2,2.2,2.4,1.3,2.3,1.7,2.7,1.5,2.5,1.6,2.6 -1 1 -2 1 $name.edgein $name.edgeout > $dir"temp"
	echo "self assignment : "`cat $dir"temp" | awk '$2 == $3' | wc -l`" out of "`cat $dir"temp" | wc -l` 
	cat $dir"temp" | awk '$2 != $3' > $name.edges
	edge_files=$edge_files$name.edges","
	rm $dir"temp"
	echo $name" total edges : "`cat $name.edges | wc -l`
	rm $dir*filtered.tsv
	
done

echo $edge_files

python -u edge_connect.py --dir $dir --edgel $edge_files --nodes $dir"ref_contiglens.tsv" --minL $minD \
--Op $Op \
--Ot $Ot \
--Om $Om \
--gff $gff

tail -n+2 $dir"summary.txt" | sed -e 's#^#'$params'#g' >> outsumm.txt 

