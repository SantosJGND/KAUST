required:
pb.fofn # list of reads fome (fasta, fastq)

modules necessary:
minimap2
samtools
python v3
	- pandas, networkx, numpy, matplotlib. 


Pipeline in order of deployment: 

## A. Map reads.

`
./purge_pipe.sh
`

change variables:

- pri_asm: contigs to bridge using reads in pb.fofn. format: fasta, extention .fa, .fasta.

Maps each **read_file** to contig file ; Processes mappings to generate table of 
- columns  : [ read, contig, length_read, length_contig, first postiion of mapping, length of mapping region, mapping quality ];
- file name: read_file.map_merge.tsv.


## B. Connect mappings - Network.

`
./net_prepare.sh
`

change variables:
- qual: [ $1 ] - minimum mapping quality [hard filter, integer] ; 
- rlen: [ $2 ] - minimum read length [hard filter, integer] ; 
- minD: [ $3 ]
- rprop: [ $4 ] - scale of contig flank size (1) [ integer ]

- Op: maximum proportion of TE [int]
- Ot: maximum number of complete TEs [int]
- Om: maximum average percentage of TEs (default 80) [int]

###
Generates output directory, named after variables above.

Filter read mappints by quality, length, and distance from contigs' starts and ends. 
Split reads into:
- .edgein: mappings to the beginning of contigs (incoming edges)
- .edgeout: mappings to the end of contigs (outgoing edges)

Merge into:
- .edges : combine incoming and outgoing edges that share same read. 

Pass all .edges files to python script *edge_connect.py*. Generates networks and identifies DAG subgraphs. 
Uses EDTA TE gff to further filter readings based on: 
	- percentage of TE overlap ; 
	- total number of entire TEs in mapped fragment (entire TE defined as 80%). 
	- average percentage of overlap with mapping (across TEs).

###

outputs to output table in current directory. 
- columns: qual, rlen, minD, rprop, totalContigs, L50, N50 (2).

(1) : relative to read length. Distance from contig start/end. 
(2) : N50 currently calculated by summing length of contigs in subgraph. This should be better, L50 is probably too large. 

