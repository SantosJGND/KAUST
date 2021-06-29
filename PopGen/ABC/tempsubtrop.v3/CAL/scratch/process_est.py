import os


import argparse

parser = argparse.ArgumentParser()

parser.add_argument("input",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read. Any number can be given.")

parser.add_argument("--temprec",type= str,default= 'filename',
	help = "slim_recipe")

args = parser.parse_args()

print()


tp_file= args.input[0]

with open(tp_file,'r') as fp:
	lines= fp.readlines()

with open(args.temprec,'r') as fp:
	temp_recipe= ''.join(fp.readlines())


lines_dict= {}

for idx in range(len(lines)):
	param= lines[idx].strip().split()

	temp_recipe= temp_recipe.replace(param[1],param[0])


new_temp_name= args.temprec.split('.')
new_temp_name= new_temp_name[:-1] + ['-temp'] + [new_temp_name[-1]]
new_temp_name= '.'.join(new_temp_name)

new_temp_name= 'instance_recipe.slim'

with open(new_temp_name,'w') as fp:
	fp.write(temp_recipe)



