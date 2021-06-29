#!/bin/bash
#SBATCH --nodes=1
#SBATCH -c 10
#SBATCH --time=120:50:00
#SBATCH --job-name=mapfq
#SBATCH --mem=80GB

module load minimap2/2.17
module load samtools/1.8


pb_list=pb.fofn
pri_asm=../../quast/Cgileadensis.hifiasm.p_ctg.2e5_new.fa_1000_0.5988_0.6430to1.8910_0.2111_primary.fasta
#pri_asm=../../haplomerge/HaploMerger2_20180603/project_CG/genome_A_ref.fa
### STEP 1. 

for i in `cat $pb_list`; 

do
	name=`basename $i` ## get read len
	sed -n '1~4s/^@/>/p;2~4p' $i > temp
	cat temp | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > $name.readlen.tsv
	sort -k1 $name.readlen.tsv > $name.readlen.sorted.tsv
	rm temp $name.readlen.tsv
	
	echo $name
	#minimap2 -t 10 -a -xmap-pb $pri_asm $i | gzip -c - > $name.sam
	samtools view -S -b $name.sam > $name.bam
	samtools sort $name.bam -o $name.sorted.bam
	samtools index $name.sorted.bam
	samtools idxstats $name.sorted.bam > $name.idxstats.txt
	sort -k1 $name.idxstats.txt > $name.idxstats.sorted.txt
	#rm $name.idxstats.txt $name.bam
	
	samtools flagstat $name.sorted.bam > $name.flagstats.out
	#samtools view $name.sam | cut -f1,3,4,5 > $name.mapIndex.tsv
	samtools view $name.sam | awk '{print $1 "\t" $3 "\t" $4 "\t" length($10) "\t" $5}' > $name.mapIndex.tsv
	sort -k2 $name.mapIndex.tsv > $name.mapIndex.sorted.tsv
	
	join -o 1.1,1.2,2.2,1.3,1.4,1.5 -1 2 -2 1 $name.mapIndex.sorted.tsv $name.idxstats.sorted.txt > $name.map_merge_int.tsv
	sort -k1 $name.map_merge_int.tsv > $name.map_merge_int.sorted.tsv	
	join -o 1.1,1.2,1.3,2.2,1.4,1.5,1.6 -1 1 -2 1 $name.map_merge_int.sorted.tsv $name.readlen.sorted.tsv > $name.map_merge.tsv
	##
	#samtools view $name.sorted.bam  | awk '{print $1"\t"length($10)}' | sort -k1 > $name.fraglen.tsv
	#sort -k1 $name.map_merge.tsv | join -o2.1,2.2,2.3,1.2,2.4,2.5,1.6,1.7 -11 -21 $name.fraglen.tsv - > $name.readmap.tsv
	
	#rm $name.sorted.bam
	echo "unmapped reads: "
	echo `samtools view -c -f4 $name.sorted.bam`
	rm $name.sorted.bam $name.bam $name.idxstats.* $name.mapIndex.tsv $name.map_merge_int.tsv 

done


