import os


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
		burnin= 4 * func_here(param_pass)

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

slim_ts = None
while slim_ts is None:
	try:
		# connect
		ts = msprime.simulate(params_ms["NANC"] * 2,
			length= params_ms["WLEN"],
			mutation_rate= mutr,
			recombination_rate= recr)
	
		slim_ts = pyslim.annotate_defaults(ts, model_type="WF", slim_generation=1)

	except:
		pass

slim_ts.dump("recipe_17.8.trees")

#####################
for param,param_pass in times_dict.items():
	temp_recipe= temp_recipe.replace(param,str(param_pass))


##################### WRITE

new_temp_name= args.temprec.split('.')
new_temp_name= new_temp_name[:-1] + ['-temp'] + [new_temp_name[-1]]
new_temp_name= '.'.join(new_temp_name)

new_temp_name= 'instance_recipe.slim'

with open(new_temp_name,'w') as fp:
        fp.write(temp_recipe)


