
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--pi', type=str, default= '',
                    help=' vcf tools windos pi file')

parser.add_argument('--mu', type=float, default= 1.3e-8,
                    help='mutation_rate')


args = parser.parse_args()

pifile=args.pi
print(pifile)
pipd= pd.read_csv(pifile, sep="\t")
picols=["CHROM","BIN_START","BIN_END","N_VARIANTS","PI"] 
print(pipd.head())

plt.figure(figsize=(10,10))
plt.hist(pipd.PI, bins= 50)

plt.ylabel('PI')
plt.title('X:{} med: {}'.format(np.mean(pipd.PI), np.median(pipd.PI)))

plt.savefig(pifile + ".pdf")

pimed=np.median(pipd.PI)

with open("NEstats.txt","a") as fp:
	net= pimed / (4 * 1.3e-8)
	nea= pimed / (4 * 6.5e-9)
	fp.write('\t'.join([str(x) for x in [pifile,pimed,net,nea]]) + "\n")
