import os
from io import StringIO
from Bio import Phylo
from Bio.Phylo.Consensus import *
import itertools as it

def read_newick(treedata):
    handle = StringIO(treedata)
    return Phylo.read(handle, "newick")

def get_handle(file):
	with open(file,'r') as fp:
		trees= fp.readlines()
	trees= [x.strip() for x in trees]
	return trees

import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--rundir', type= str, default= "run1/")
parser.add_argument('-o', type= str, default= "majority_tree.nex")
parser.add_argument('-m', type= str, default= "newick")

args = parser.parse_args()

print(args.rundir)
rundir=args.rundir

content= os.listdir(rundir)

trees= [x for x in content if "treeout" in x]
trees= [get_handle(rundir + x) for x in trees]
trees= list(it.chain(*trees))
trees= [read_newick(x) for x in trees]

majority_tree= majority_consensus(trees, 0.5)
support_tree = get_support(majority_tree, trees)

Phylo.draw_ascii(support_tree)
Phylo.write(support_tree, args.o, args.m)
