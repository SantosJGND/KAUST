
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--taj', type=str, default= '',
                    help=' vcf tools windos pi file')

args = parser.parse_args()

pifile=args.taj
print(pifile)
pipd= pd.read_csv(pifile, sep="\t")
picols=["CHROM","BIN_START","N_SNPS","TajimaD"] 
print(pipd.head())

plt.figure(figsize=(10,10))
plt.hist(pipd.TajimaD, bins= 50)

plt.ylabel('TajimaD')
plt.title('X:{} med: {}'.format(np.mean(pipd.TajimaD), np.nanmedian(pipd.TajimaD)))

plt.savefig(pifile + ".pdf")
