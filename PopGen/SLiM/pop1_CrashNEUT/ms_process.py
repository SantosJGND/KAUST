import numpy as np
import os
import collections

def write_ms(
    tree_sequence,
    output,
    print_trees=False,
    precision=4,
    num_replicates=1,
    write_header=True,
):
    """
    Write ``ms`` formatted output from the genotypes of a tree sequence
    or an iterator over tree sequences. Usage:
    .. code-block:: python
        import tskit as ts
        tree_sequence = msprime.simulate(
            sample_size=sample_size,
            Ne=Ne,
            length=length,
            mutation_rate=mutation_rate,
            recombination_rate=recombination_rate,
            random_seed=random_seed,
            num_replicates=num_replicates,
        )
        with open("output.ms", "w") as ms_file:
            ts.write_ms(tree_sequence, ms_file)
    :param ts tree_sequence: The tree sequence (or iterator over tree sequences) to
        write to ms file
    :param io.IOBase output: The file-like object to write the ms-style output
    :param bool print_trees: Boolean parameter to write out newick format trees
        to output [optional]
    :param int precision: Numerical precision with which to write the ms
        output [optional]
    :param bool write_header: Boolean parameter to write out the header. [optional]
    :param int num_replicates: Number of replicates simulated [required if
        num_replicates used in simulation]
    The first line of this ms-style output file written has two arguments which
    are sample size and number of replicates. The second line has a 0 as a substitute
    for the random seed.
    """
    if not isinstance(tree_sequence, collections.abc.Iterable):
        tree_sequence = [tree_sequence]

    i = 0
    for tree_seq in tree_sequence:
        if i > 0:
            write_header = False
        i = i + 1

        if write_header is True:
            print(
                f"ms {tree_seq.sample_size} {num_replicates}",
                file=output,
            )
            print("0", file=output)

        print(file=output)
        print("//", file=output)
        if print_trees is True:
            """
            Print out the trees in ms-format from the specified tree sequence.
            """
            if len(tree_seq.trees()) == 1:
                tree = next(tree_seq.trees())
                newick = tree.newick(precision=precision)
                print(newick, file=output)
            else:
                for tree in tree_seq.trees():
                    newick = tree.newick(precision=precision)
                    print(f"[{tree.span:.{precision}f}]", newick, file=output)

        else:
            s = tree_seq.get_num_sites()
            print("segsites:", s, file=output)
            if s != 0:
                print("positions: ", end="", file=output)
                positions = [
                    variant.position / (tree_seq.sequence_length)
                    for variant in tree_seq.variants()
                ]
                for position in positions:
                    print(
                        f"{position:.{precision}f}",
                        end=" ",
                        file=output,
                    )
                print(file=output)

                genotypes = tree_seq.genotype_matrix()
                print(genotypes.shape)
                print(tree_seq.num_samples)
                gensamp= np.array(genotypes,dtype= str).T
                print("\n".join(["".join(x) for x in gensamp]), file= output)
                #for k in range(tree_seq.num_samples):
                #    tmp_str = "".join(map(str, genotypes[:, k]))
                #    print(tmp_str, file=output)

            else:
                print(file=output)


import argparse

parser = argparse.ArgumentParser()

parser.add_argument("input",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read. Any number can be given.")

parser.add_argument("--temprec",type= str,default= 'filename',
        help = "slim_recipe")

parser.add_argument("--rescale",type= float, default= 1.0, 
	help= 'rescale parameters if action != none.')

parser.add_argument("--burnin",type= int,default= 4e4,
        help ="burnin in gens.")

parser.add_argument("--out",type=str,default= "sims/",
        help ="out_directory")

parser.add_argument("--tag",type= str,default= "0",
        help ="tag")


args = parser.parse_args()

burnin= args.burnin

tp_file= args.input[0]

with open(tp_file,'r') as fp:
        lines= fp.readlines()

with open(args.temprec,'r') as fp:
        temp_recipe= ''.join(fp.readlines())

times_dict= {}
lines_dict= {}

func_operate= [float,int]

params_ms= { x:0 for x in ["NANC","WLEN","RECR","MUTR"] }


for idx in range(len(lines)):
	if len(lines[idx]) < 2:
		continue
	param= lines[idx].strip().split()
	print(param)

	func_here= func_operate[int(param[2])]

	action= ""
	if len(param) == 4:
		action= param[3]

		if action== "+":
			param_pass= float(param[0]) * args.rescale

		if action=="r":
			param_pass= .5 * (1 - (1 - 2*float(param[0]))**args.rescale)

		if action=="none":
			param_pass= param[0]

		if action== "b":
			param_pass = float(param[0]) / args.rescale
	else:
		param_pass= float(param[0]) / args.rescale

	param_pass= func_here(param_pass)

	if action != 'b':
		temp_recipe= temp_recipe.replace(param[1],str(param_pass))
	else:
		times_dict[param[1]]= param_pass
	
	if param[1] == "NANC":
		NANC= func_here(param_pass)
		burnin= 10 * func_here(param_pass)

	if param[1] in params_ms.keys():
		params_ms[param[1]]= func_here(param_pass)
	
	if param[1] == "SRW":
		params_ms[param[1]]= func_here(param_pass)


######
###### BURN-IN MS 

import msprime, pyslim

print(params_ms["NANC"])
print(params_ms["WLEN"])
print(params_ms["MUTR"])
print(params_ms["RECR"])

mutr=params_ms["MUTR"]
recr=params_ms["RECR"]

if "SRW" in params_ms:
	F= params_ms["SRW"] / (2 - params_ms["SRW"])
	
	mutr= mutr / (1 + F)
	recr= recr * (1 - params_ms["SRW"])

print('ms rec: {}'.format(recr))
print('ms mut: {}'.format(mutr))

slim_ts = None
#while slim_ts is None:
#	try:
		# connect

print("#")
print(params_ms["NANC"])
ps = msprime.simulate(
	Ne= params_ms["NANC"],
	sample_size= params_ms["NANC"] * 2,
	length= params_ms["WLEN"],
	mutation_rate= mutr,
	recombination_rate= recr,
	num_replicates= 1)

#slim_ts = pyslim.annotate_defaults(ps, model_type="nonWF", slim_generation=1)


print("hi")
tree_dump= "ANC_trees/" + args.tag + ".trees"
ms_dump= "ANC_trees/" + args.tag + ".ms"
#slim_ts.dump(tree_dump)

import tskit as ts
#ts1=ts.load(tree_dump)

for tree in ps:
	with open(ms_dump,"w") as ms_dump:
    		write_ms(tree, ms_dump, precision=6, num_replicates=1, write_header=True)

temp_recipe= temp_recipe.replace("msdump", '"{}"'.format(ms_dump))
#####################
temp_recipe= temp_recipe.replace("BURNIN",str(burnin))
for param,param_pass in times_dict.items():
	temp_recipe= temp_recipe.replace(param,str(param_pass + burnin))


##################### WRITE
outms= args.out + args.tag + ".ms"
outfixed= args.out + args.tag + ".fixed"

temp_recipe= temp_recipe.replace("outfilems", '"{}"'.format(outms))
temp_recipe= temp_recipe.replace("outfixed", '"{}"'.format(outfixed))
temp_recipe=temp_recipe.replace("treefile", '"{}"'.format(tree_dump))

####################

new_temp_name= args.temprec.split('.')
new_temp_name= new_temp_name[:-1] + ['-temp'] + [new_temp_name[-1]]
new_temp_name= '.'.join(new_temp_name)

new_temp_name= args.out + args.tag + '.slim'

with open(new_temp_name,'w') as fp:
        fp.write(temp_recipe)


