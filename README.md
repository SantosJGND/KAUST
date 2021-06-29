## sub-projects scripts

### bridge_contigs/

Map reads to assembly, find reads that bridge contigs.Find *contig connected components*; summary statistics.

work directory:
/ibex/scratch/projects/c2060/SANTOSJ_DIR/Projects/Commiphora/assembly_double/assembly_refinements/mappings

### Distribution/

Kolmogorov-Smirnof test across contigs, readme inside. 

### GeneXtract_v1/

Map and fasta to list of assemblies, extract best positions and corresponding sequences, 
- subset orfs; 
- predict protein;
- construct graph; 

### IIMs/

Calculate ld within inversions, use Inversion matrix as input. 

- working directory: 
/ibex/scratch/projects/c2060/SANTOSJ_DIR/subProjects/IIMs


### merge_PAVs/

Merge presence absence calls from different software. Given standardized input format. 

i. merge Insertions and Deletions separately.
ii. Produce overlap statistics and 3-way venn diagrams.

### Ocoarcata_demography/

Demographic inference using single nucleotide polimorphism data and fastsimcoal2.6 software.

- Three population scenarios, with and without migration;
- Setup for independent deployment across subgenomes.

Used for inference of Oryza Coarctata.

- data directory:
/ibex/scratch/projects/c2060/SANTOSJ_DIR/Projects/Ocoarc/

### PopGen/

Population genetics applications. 

- Fst calculatons, 
- Window-based estimates of diversity (pi and tajD), 
- Window-based conversion from plink tyo; 
- treemix replicates consensus tree extraction ; 
- Demogaphic simulations using SLiM - alternative 2 population scenarios with and without positive and background selection. 
- ePSMC demographic inference and selfing-rate estimation ; 
- ABC inference;

### INVtype/

This was meant to determine inversion orientatio from mapping strands. But it was a bust. Still.

Given region in reference assembly (samtools format), and mapping bamfile (against reference assembly):
- extract alignments to selecte region.
- re-map against extract fasta to generate paf.

- readme.md inside. 


