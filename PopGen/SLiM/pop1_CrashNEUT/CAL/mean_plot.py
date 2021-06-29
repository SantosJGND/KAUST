import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



file= "simCrashSel_200.stats"

dbstat_source= pd.read_csv(file,sep= "\t",error_bad_lines=False, warn_bad_lines= False)

param_dict= {
        "NTP0": [10000,20000,40000],
        "SRD": [.99],
        "NCD": [600]
#        "P0": [0.001,0.5],
#        "Fs": [0.01,0.5],
}

compare= [x for x,g in param_dict.items() if len(g) >1]

stats= ["thetapi", "thetaw","thetah", "hprime","tajimasd","hapdiv"]

stat_collect= [0,1,2,4,5]

stat_borders= {
	"hprime": [-0.05,1]
        "thetaw": [-.05,1],
        "thetapi": [-.05,1],
        "tajimasd": [-.8,.8],
        "thetah": [-.05,1],
        "hapdiv": [0,1]
}


###
total_length= 52560

###
dict_obs= {}

import itertools as it
param_keys= list(param_dict)
comp_idx= [param_keys.index(x) for x in compare]

param_values= [param_dict[g] for g in param_keys]
param_combs= list(it.product(*param_values))

for comb in param_combs:
	dbstat= dbstat_source.copy()

	for idx in range(len(comb)):

		param= param_keys[idx]
		g= comb[idx]
		print(param, g)
		dbstat= dbstat[dbstat[param] == g]
        
	comb_str= "-".join([str(x) for x in comb])
	dict_obs[comb_str]= dbstat

	## merge windows
	wind_obs= {}
	
	for wind in dbstat.win.unique():
		wind_sel= dbstat[dbstat.win == wind]

		wind_obs[wind]= {
			"mean": wind_sel.mean(axis= 0),
			"sd": wind_sel.mean(axis= 0)
		}
	dict_obs[comb_str]= wind_obs

fixed= ["".join([v,str(g[0])]) for v,g in param_dict.items() if len(g) == 1]
fixed= ".".join(fixed)

for stat in stat_borders.keys():
	
	plt.figure(figsize= (10,10))
	
	data= []
	data_names= []
	for comb,g in dict_obs.items():
		combstats= comb.split("-")
		stat_part= [param_keys[x] + ":" + combstats[x] for x in comp_idx]
		gstat= [dict_obs[comb][w]["mean"][stat] for w in dict_obs[comb].keys()]
		data.append(gstat)
		data_names.append(stat_part)
		#plt.hist(g[stat],label= comb, density= True, alpha= .6)
	
	plt.boxplot(data, labels= data_names)
	plt.title(stat)
	plt.ylim(tuple(stat_borders[stat]))
	fig_name= fixed + "_" + "-".join(compare) + "_" + stat + ".pdf"
	plt.savefig(fig_name)
