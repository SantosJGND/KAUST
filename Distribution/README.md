# Uniform distribution test


Given a bed file, the script Distribution_simple.py performs Kolmogorov-Smirnov (KS) tests of uniform distribution per contig. 
The False positive rate of this test is estimated by performing Ni simulations of randomly simulated uniform distributions, given the same sample size. 

The point of this control is to study the impact of sample size. 

> KS tests are performed on the middle coordinate of each segment (TE). 

## running.

The script takes the following arguments:
--input : a bed file ; 
--alpha : significance threshold, to estimate significance of individual tests. 

the arguments --ksteps and --ksteps are used to determine the false positive rate (kwindow is no longer necessary). 

## output:
table in directory [input_alpha]; 
columns: CHROM   NINV    KS      PVAL    FDR

FDR should be ~alpha ; 

The PVAL column should be divided by the number of tests performed (the number of contigs, i.e. Bonferroni correction). 

