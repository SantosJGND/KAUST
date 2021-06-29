import os
from os import path
import numpy as np

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("input",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read. Any number can be given.")

args = parser.parse_args()



vcf_file= args.input[0]

if len(args.input) == 1:
	out_file= "summary_stats_temp.txt"
else:
	out_file= args.input[1]

with open(vcf_file,'r') as fp:
	lines= fp.readlines()
	lines= [x.strip().split() for x in lines]

lines= lines[2:]
lines= np.array(lines,dtype= str)
lines= lines.T[1:]

with open(out_file,'w') as fp:
	for i in range(lines.shape[0]):
		for j in range(lines.shape[1]):
			fp.write('f{}_{}'.format(i,j) + '\t')
	fp.write('\n')
	for i in range(lines.shape[0]):
		for j in range(lines.shape[1]):
			fp.write(str(lines[i,j]) + '\t')

