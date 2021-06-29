
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt



vcf_file="out.weir.fst"

fst_array= pd.read_csv(vcf_file,sep= "\t")
fst_array=fst_array.dropna()

fst_array=fst_array.loc[fst_array.WEIR_AND_COCKERHAM_FST > -.2]

figdir= "Figures/"

plt.figure(figsize= (8,8))

X= list(fst_array["WEIR_AND_COCKERHAM_FST"])

print(np.mean(X))

plt.hist(X,20,density=True, alpha=0.75)

plt.title("mean: {}".format(np.mean(X)))

plt.xlim(-.2,1)
plt.ylabel('density')
plt.xlabel('Weir-Cockerham FST')

plt.savefig(figdir + 'hist.pdf')
