## Genome Size Estimation using GenomeScope 

### installing GenomeScope;


1. load R, v3.6.0

2. do `echo $R_LIBS` to know in what directory your R libraries should go in. Create it if necessary. 

3. Inside **install.R**, replace *local_lib_path* by your $R_LIBS. 

4. run `Rscript install.R`


> Problems: how i solved my minpack.lm installation:

i. remove the directory:
	- $R_LIBS/00LOCK-minpack.lm


### Kmer Histogram

1. Modify the **JellyGS_demo.sh**

check kmer size and input (directory + filename, note the wildcard used to represent more than one file).

`
./jellyGS_demo.sh
`


### GenomeScope

1. modify the Script **GS_demo.sh**; check that kmer size is the same as used to produce the kmer histogram. 

`
./GS_demo.sh
`

## Exercises:

i. Collapse the two scripts, **jellyGS_demo.sh** and **GS_demo.sh** into a single one. 
	- pay attention to modules that need loading./ 


ii. Make kmer size into and variable passed to the script from the command line (i.e. kmer=$1).

	- vary kmer size. Study output. 
