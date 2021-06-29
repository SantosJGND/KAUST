import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



file= "simCrashSel_200.stats"

dbstat_source= pd.read_csv(file,sep= "\t",error_bad_lines=False, warn_bad_lines= False)

param_dict= {
        "NTP0": [20000],
        "SRD": [0,.8],
        "NCD": [600,5000],
        "P0": [0.001,0.5],
        "Fs": [0.01,0.5],
}

stats= ["thetapi", "thetaw","thetah", "hprime","tajimasd","hapdiv"]

stat_collect= [0,1,2,4,5]

stat_borders= {
	"thetaw": [-.05,.8],
        "thetapi": [-.05,.8],
        "tajimasd": [-.8,.8],
        "thetah": [-.05,.8],
	"hapdiv": [0,.8]
}


###
total_length= 52560

###
import itertools as it
param_keys= list(param_dict)
param_values= [param_dict[g] for g in param_keys]
param_combs= list(it.product(*param_values))

for comb in param_combs:
	dbstat= dbstat_source.copy()

	for idx in range(len(comb)):

		param= param_keys[idx]
		g= comb[idx]
		print(param, g)
		dbstat= dbstat[dbstat[param] == g]


	print(dbstat.head())

	## merge windows
	dict_obs= {}

	for wind in dbstat.win.unique():
		wind_sel= dbstat[dbstat.win == wind]

		dict_obs[wind]= {
			"mean": wind_sel.mean(axis= 0),
			"sd": wind_sel.mean(axis= 0)
		}


	###
	###
	for idx in stat_collect:

		stat= stats[idx]
		print(stat)

		pltname= "-".join([".".join([param_keys[v],str(comb[v])]) for v in range(len(comb))])
		pltitle= "; ".join([": ".join([param_keys[v],str(comb[v])]) for v in range(len(comb))])

		plt.figure(figsize= (15,10))

		X= sorted(dict_obs.keys())
		Y= [dict_obs[x]["mean"][stat] for x in X]

		yerr= [dict_obs[x]["sd"][stat] for x in X]
		X= np.array(X) * (total_length / len(dict_obs))

		print("#")

		plt.plot(X,Y)
		plt.vlines(50000,stat_borders[stat][0],stat_borders[stat][1])

		plt.xlabel("windows", fontsize= 20)
		plt.ylabel(stat, fontsize= 20)
		plt.ylim(stat_borders[stat][0],stat_borders[stat][1])
		plt.title(pltitle + ' mean {} : {}'.format(stat,np.mean(Y)), fontsize= 15)

		plt.savefig(pltname + "_{}.pdf".format(stat))
		plt.close()







