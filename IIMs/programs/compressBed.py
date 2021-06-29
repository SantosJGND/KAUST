import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--file', type=str,
                    default='file.tsv')


parser.add_argument('--out', type=str,
                    default='')

args = parser.parse_args()


outfile= args.out

if not args.out:
        outfile= args.file


with open(args.file,'r') as fp:
        lines= fp.readlines()

lines= [x.split() for x in lines]

chroms= [x[0] for x in lines]
chroms_idx= {z: [x for x in range(len(lines)) if chroms[x] == z] for z in list(set(chroms))}

chroms_dict= {
        z: {
                's': [int(lines[x][1]) for x in g],
                'e': [int(lines[x][2]) for x in g]
        } for z,g in chroms_idx.items()
}

with open(outfile,'w') as fp:
        for chrom,g in chroms_dict.items():
                nline= '{}:{}-{}'.format(chrom,min(g['s']),max(g['e']))
                fp.write(nline + '\n')


