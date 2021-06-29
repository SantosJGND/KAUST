import Bio.SeqIO as IO
import gzip
import pandas as pd

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("--refseq", type= str, default= "GCF_001433935.1_IRGSP-1.0_genomic.fna.gz",
			help="assembly fasta, gzipped")


parser.add_argument("--reg", type= str, default= "Chr01:0-10000",
                        help="fasta region to get")

parser.add_argument("-o", type= str, default= "output.fa",
                        help="output file")

args = parser.parse_args()

### Process region
###
chrom_file="/ibex/scratch/santj0a/GenomeGraphs/VG/GENE_collect/prelim_data/IRGSP_chroms.txt"

chrom, interval= args.reg.split(":")
start, end= [int(x) + 1 for x in interval.split("-")]

if "IRGSP" in args.refseq:
	chromdb=pd.read_csv(chrom_file,sep= "\t")
	chromsect=chromdb[chromdb.default == chrom].irgsp[0]
	print("chainging chrom {} to {}".format(chrom, chromsect))
	chrom= chromsect
###
### Get sequence
seqA_file= args.refseq

with gzip.open(seqA_file,"rt") as fp:
        seqA= IO.to_dict(IO.parse(fp, "fasta"))

regA= seqA[chrom]
seqA= str(regA.seq[start:end])

###
### Write file

with open(args.o, 'w') as fp:
	fp.write('>{} {}-{}'.format(chrom, start, end) + "\n")
	fp.write(seqA)
	

