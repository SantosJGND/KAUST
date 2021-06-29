## Inversions Informative Markers

Given an Inversion matrix (first four columns = bed fle), and plink file; 

calculate linkage disequilibrium among accessions marked as Inverted in inversion matrix. 

The python sript `INV_select.py` parses matrix to determine which accessions are marked as Inverted for a particular selected inversion. 

plink is used to subset data (missingness, maf) and calculate ld. 

programs/parse_snps.py is used to parse the final output for a differences in ld and threshold frequency.

### deployment

i. Prepare Inversion Matrix in IIMs/prelim_data/INV_matri.txt.

ii. Modify `programs/parse_snps.py` for ld and frequency thresholds.

iii. With index = the index of the wanted inversion in matrix; id= run identifier.

`
./INV_dispatch.sh $index $id

`

