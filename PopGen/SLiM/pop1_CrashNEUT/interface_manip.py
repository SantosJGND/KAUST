import numpy as np

interface= "AB-SinterfaceR.txt"

def_names= ["$f0","$f1","$f2","$f3"]
params_dfe= {
	"DF1": [0.7,0.1,0.1,0.1],
	"DF2": [0.1,0.1,0.1,0.7],
	"DF3": [0.25,0.25,0.25,0.25]
}

params_size= {
	"$NANC": [100000],
	"$NTP0": [10000, 20000, 40000],
	"$SRW": [0.55],
	"$SRD": [0, 0.5, 0.8, 0.99],
	"$NCD": [600, 5000]
}


with open(interface,'r') as fp:
	interstore= "".join(fp.readlines())

import itertools as it
home= "ABR_params/"

combs= [params_size["$NANC"],params_size["$NTP0"],params_size["$SRW"],params_size["$SRD"],params_size["$NCD"]]

list_combs= list(it.product(*combs))

for param_comb in list_combs:
	new_interface= str(interstore)
	param_comb= [str(x) for x in param_comb]
	new_interface= new_interface.replace("$NANC", param_comb[0])
	new_interface= new_interface.replace("$NTP0", param_comb[1])
	new_interface= new_interface.replace("$SRW", param_comb[2])
	new_interface= new_interface.replace("$SRD", param_comb[3])
	new_interface= new_interface.replace("$NCD", param_comb[4])
	
	combname= '-'.join(np.array(param_comb,dtype= str)) + '_interface.txt'
	
	with open(home + combname,'w') as fp:
		fp.write(new_interface)
	
