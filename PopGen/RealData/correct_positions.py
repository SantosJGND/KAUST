import pandas as pd
import numpy as np

import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--ms', type=str,
                    help='an integer for the accumulator')

parser.add_argument('--windows', type=str,
                    help='an integer for the accumulator')

parser.add_argument('--wlen', type=int, default= 100000,
                    help='length of region')


args= parser.parse_args()

ms_file= args.ms
window_bed=args.windows

with open(ms_file,'r') as fp:
	lines= fp.readlines()
	
	where_pos= [x for x in range(len(lines)) if 'positions' in lines[x]][0]
	pos= lines[where_pos]
	pos= pos.strip("\n").strip("positions:").split()

#print(pos)

pos= np.array(pos,dtype= int)

band= pd.read_csv(window_bed, sep= "\t", header= None)
band.columns= ["chrom","IN","OUT","ID"]
band["height"] = band.OUT - band.IN
total= np.sum(band.height)
print("height: {}".format(total))

#print(band)

####
####
start= band.IN[0]
spaces= [0]
for idx in range(1,band.shape[0]):
	gulf= band.IN[idx] - band.OUT[idx-1]
	spaces.append(gulf)

#print("spaces: {}".format(spaces))

new_pos= []
for idx in range(len(pos)):
	posi= pos[idx]
	where= [x for x in range(band.shape[0]) if band.OUT[x] >= posi][0]
	new_pos.append(posi - start - sum(spaces[:(where+1)]))

new_pos= np.array(new_pos,dtype= float) / total
#print(new_pos)
new_pos= np.array(new_pos,dtype= str)

new_line= "positions: " + " ".join(new_pos) + "\n"
lines[where_pos]= new_line

with open(ms_file[:-3] + "_compressed.ms","w") as fp:
	fp.write("".join(lines))

