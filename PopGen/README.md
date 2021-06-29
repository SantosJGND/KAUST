### PopGen/

Population genetics applications.


### FST/
- Fst calculatons,

scripts:
> fst_calc.sh; plot_launch.sh;

### PopGenStats/

- Window-based estimates of diversity (pi and tajD) from VCF file.
	- takes repeat mask file. 
	- calculations per sub-group: see IDs/ directory.
scripts:
> pi_get.sh

### ePSMC/


### treemix/

Extracting consensus tree from treemix replicates in python. 

### RealData/

Given plink file and list of subset IDs (underscore delimited); 
	- extract local genomic windows of given size, convert to MS format. 

*used to compare with SLiM simulations in SLiM/; 

### SLiM/

Demogaphic simulations using SLiM - alternative 2 population scenarios with and without positive and background selection.

subdirectories:
- `CAL/`: calculate statistics per simulation ; 
- `ePSMC`: demographic inference on simulated dta. 
- `recipes`: slim recipes: demographic scenarios and genome structure (gene density)
- `environments`: yaml files for environments used.
- `pop1_CrashNEUT`: example run directory. 

### ABC/

Approximate Bayesian Computation: determining best demographic parameters for single model for future model compaison. 

subdirectories: 
- `programs`: copy of main scripts, utilities scripts. 
- `tempsubtrop.v3` : example run. 


