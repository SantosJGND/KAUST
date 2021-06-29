
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import itertools as it

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)


class db_stats():
        def __init__(self, db_source):
                self.db= pd.read_csv(db_source,sep= "\t",error_bad_lines=False, warn_bad_lines= False)

        def params_find(self,param_dict):
                self.param_dict= param_dict
                self.keys= list(param_dict.keys())
                fixed= ["".join([v,str(g[0])]) for v,g in param_dict.items() if len(g) == 1]

                self.fixed= ".".join(fixed)
                self.var= [x for x,g in param_dict.items() if len(g) >1][0]
                self.varidx= self.keys.index(self.var)

        def summarize(self, stats):
                dict_obs= {}

                param_keys= self.keys
                param_values= [param_dict[g] for g in param_keys]
                param_combs= list(it.product(*param_values))

                for comb in param_combs:
                        dbstat= self.db.copy()
                        stat_part= self.var + ":" + str(comb[self.varidx])

                        for idx in range(len(comb)):
                                param= param_keys[idx]
                                g= comb[idx]
                                #print(param, g)
                                dbstat= dbstat[dbstat[param] == g]

                        comb_str= "-".join([str(x) for x in comb])
                        #dict_obs[comb_str]= dbstat

                        ## merge windows
                        wind_obs= {}

                        for wind in dbstat.win.unique():
                                wind_sel= dbstat[dbstat.win == wind]

                                wind_obs[wind]= {
                                        "mean": wind_sel.mean(axis= 0),
                                        "sd": wind_sel.mean(axis= 0)
                                }

                        gstats= {
                                stat: [wind_obs[x]["mean"][stat] for x in wind_obs.keys()] for stat in stats
                        }

                        print(gstats.keys())

                        dict_obs[stat_part]= gstats

                self.chain= dict_obs
                self.stats= stats

        def boxplot(self, stat_borders= {}):

                for stat in self.stats:

                        plt.figure(figsize= (10,10))

                        data= []
                        data_names= []

                        for stat_part,gstat in self.chain.items():
                                data.append(gstat[stat])
                                data_names.append(stat_part)

                        plt.boxplot(data, labels= data_names)
                        plt.title(stat)

                        if stat in stat_borders.keys():
                                plt.ylim(tuple(stat_borders[stat]))

                        fig_name= self.fixed + "_" + self.var + "_" + stat + ".pdf"
                        plt.savefig(fig_name)

#####################
#####################
#####################

stats= ["thetapi", "thetaw","thetah", "hprime","tajimasd","hapdiv"]

stat_borders= {
        "hprime": [-1,.5],
        "thetaw": [-.05,.5],
        "thetapi": [-.05,.6],
        "tajimasd": [-1,1],
        "thetah": [-.05,1],
        "hapdiv": [0,.4]
}

################################################################################
################################################################################

file= "/ibex/scratch/santj0a/Projects/SLiM/RealData/CAL/simCrashSel_200.stats"

dbstat_source= pd.read_csv(file,sep= "\t",error_bad_lines=False, warn_bad_lines= False)

###
total_length= 100000
###
dbstat= dbstat_source.copy()

## merge windows
wind_obs= {}

for wind in dbstat.win.unique():
        wind_sel= dbstat[dbstat.win == wind]

        wind_obs[wind]= {
                "mean": wind_sel.mean(axis= 0),
                "sd": wind_sel.mean(axis= 0)
        }

dict_obs= {
        stat: [wind_obs[x]["mean"][stat] for x in wind_obs.keys()] for stat in stats
}


##################################
##################################
title= "NEUTSRD"

file= "/ibex/scratch/santj0a/Projects/SLiM/pop1_NoCrashNEUT/CAL/simCrashSel_200.stats"

param_dict= {
        "NTP0": [20000],
        "SRW": [0.55],
        "SRD": [0,.8,.99]
}

crash_demo= db_stats(file)
crash_demo.params_find(param_dict)
crash_demo.summarize(stats)
#crash_demo.boxplot()


##################################
##################################
file= "simCrashSel_200.stats"

param_dict= {
        "NTP0": [20000],
        "SRD": [0,.8,.99],
        "NCD": [600]
}

small_crash_demo= db_stats(file)
small_crash_demo.params_find(param_dict)
small_crash_demo.summarize(stats)
#crash_demo.boxplot()


##################################
##################################
label0= "real_data"
label1= crash_demo.fixed
label2= small_crash_demo.fixed

for stat in stats:

        ticks= list(crash_demo.chain.keys())
        data_rd= [dict_obs[stat]] * len(ticks)
        data_a= [crash_demo.chain[x][stat] for x in ticks]
        data_b= [small_crash_demo.chain[x][stat] for x in ticks]

        plt.figure()
        bp0 = plt.boxplot(data_rd, positions=np.array(range(len(data_rd)))*3.0-0.7, sym='', widths=0.6)
        bpl = plt.boxplot(data_a, positions=np.array(range(len(data_a)))*3.0, sym='', widths=0.6)
        bpr = plt.boxplot(data_b, positions=np.array(range(len(data_b)))*3.0+0.7, sym='', widths=0.6)
        set_box_color(bp0, '#000')
        set_box_color(bpl, '#D7191C') # colors are from http://colorbrewer2.org/
        set_box_color(bpr, '#2C7BB6')

        # draw temporary red and blue lines and use them to create a legend
        plt.plot([], c='#000', label=label0)
        plt.plot([], c='#D7191C', label=label1)
        plt.plot([], c='#2C7BB6', label=label2)
        plt.legend()
        plt.title(stat,fontsize= 15)
        for vl in [1.5, 4.5]:
            plt.axvline(vl, color = 'r', linestyle='--')

        plt.xticks(range(0, len(ticks) * 3, 3), ticks)
        plt.xlim(-1.6, len(ticks)*(len(ticks)-1) + 1.6)
        plt.ylim(tuple(stat_borders[stat]))
        plt.tight_layout()
        plt.savefig('{}_compare_{}.png'.format(stat,title))


