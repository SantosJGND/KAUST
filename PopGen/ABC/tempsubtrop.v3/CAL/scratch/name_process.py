import argparse

parser = argparse.ArgumentParser()

parser.add_argument("input",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read. Any number can be given.")


parser.add_argument("--tag",type= str,default= 'filename',
	help = "new col")

parser.add_argument("--sep",type= str,default= '\t',
        help = "sep")

args = parser.parse_args()

for file in args.input:
	with open(file,'r') as fp:
		lines= fp.readlines()
	
	lines= [x.strip() + args.sep + args.tag for x in lines]
	
	with open(file,'w') as fp:
		fp.write('\n'.join(lines) + '\n')

