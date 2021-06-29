#!/bin/bash

outf=summary.txt

echo -e "folder\tchrom\tstart\tend\tlen\tNref\tNcontrol\tIIMs" > $outf

for dir in run*/; do 

if [ -f $dir"base_info/region.txt" ]; then 

region=`cut -f1-3 $dir"base_info/region.txt"`
region_len=`cat $dir"base_info/region_length.txt"`
controln=`cat $dir"base_info/control_fasta.txt" | wc -l`
refn=`cat $dir"base_info/ref_fasta.txt" | wc -l`

folder=`basename $dir`
nline=$folder"\t"$region"\t"$region_len"\t"$refn"\t"$controln

if [ -f $dir"selected.snps" ]; then 

nsel=`cat $dir"selected.snps" | wc -l`
nline=$nline"\t"$nsel

else
nline=$nline"\tNA"
fi

echo -e $nline >> $outf

fi
done

sed -i -e 's#[[:space:]]#\t#g' $outf
