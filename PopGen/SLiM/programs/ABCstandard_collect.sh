#!/bin/bash

range=20
file_in=tempsubtrop_sim_sampling1.txt
file_out=tempsubtrop_sim_obs0_sampling1.txt

d=0

for i in `seq 1 $range`; do
	nfile=scratch$i"/"$file_in
	if [ $d = 0 ]; then
		cat $nfile > $file_out
		d=1
	else
		tail -n +2 $nfile >> $file_out
fi

done
