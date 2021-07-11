## Merge gff

Input two gff files, gff1 & gff2. 

i. filter for tag in third column
	- default: "gene" for gff1; LTR for gff2 ; 

assumes format:
	["contig","m1","pred","start","end","q","s","t","idt"]
	- m1, q,s,t are placeholders. columns not used. important is 9 columns and position of contig, pred, start, end and idt. 

ii. overlap and distance.

for every element in gff1:
	- extract overlaping elements in gff2;
	- extract distance to closes element in gff2;

iii. Output: 

summary table.
Columns:
	- contig ;
	- start position ; 
	- end position ;
	- id (gff1) 
	- id of closest element (gff2)
	- distance to closest element (gff2, base pairs) ; 
	- overlap - number of overlapping gff2 elements ;
	- overlap_id - id list of overlapping gff2 elements ; 
