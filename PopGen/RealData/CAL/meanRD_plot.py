import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



file= "/ibex/scratch/santj0a/Projects/SLiM/RealData/CAL/simCrashSel_200.stats"

dbstat_source= pd.read_csv(file,sep= "\t",error_bad_lines=False, warn_bad_lines= False)

param_dict= {
        "scratch": []
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
total_length= 100000

###

dbstat= dbstat_source.copy()


dict_obs= {}

## merge windows
wind_obs= {}

for wind in dbstat.win.unique():
        wind_sel= dbstat[dbstat.win == wind]

        wind_obs[wind]= {
                "mean": wind_sel.mean(axis= 0),
                "sd": wind_sel.mean(axis= 0)
        }

dict_obs["RD"]= wind_obs

plt.figure(figsize= (10,10))

comb= "RD"
g= dict_obs[comb]

data= []
data_names= []

for stat in stat_borders.keys():

        combstats= comb.split("-")
        stat_part= stat
        gstat= [dict_obs[comb][w]["mean"][stat] for w in dict_obs[comb].keys()]
        data.append(gstat)
        data_names.append(stat_part)
        #plt.hist(g[stat],label= comb, density= True, alpha= .6)

plt.boxplot(data, labels= data_names)
plt.title(stat)
fig_name= "RDstats.pdf"
plt.savefig(fig_name)


