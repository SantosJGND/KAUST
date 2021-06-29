import pandas as pd

summary= "summary.txt"
threshold= 0.4
min_diff= 0.3

summ= pd.read_csv(summary,sep= "\t")
print(summ.shape)

summ["diffs"]= abs(summ.ALT1 - summ.ALT2)
print(summ.head())

summ= summ[(summ.ALT1 > 0) | (summ.ALT2 > 0)].reset_index(drop= True)
summ= summ[(summ.diffs > min_diff)].reset_index(drop= True)
summ= summ[summ.WEIR_AND_COCKERHAM_FST > threshold].reset_index(drop=True)
print(summ.head())

snps= list(summ.SNP)
snps= [str(w) for w in snps]

summ.to_csv("selected.snps", sep= "\t")

with open("snplist",'w') as fp:
	fp.write("\n".join(snps))
