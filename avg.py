from multiprocessing import Pool
import numpy as np
import os, time
import sys
import pickle as pkl


N_avo=6.0221409e23
vol_er = 3.9*0.1*0.1

#get data location
with open("outputLoc.mdl",'r') as f:
	exec(f)

#set mcell directory
dir=outputLoc+"ppf/"+sys.argv[1]+"/"

#averaged data folder
avg_path=os.path.join(dir,"Average")
if not os.path.exists(avg_path):
    os.makedirs(avg_path)
    
#seed folders
seed_folders=os.listdir(dir)
seed_folders.remove("Average")

#list of file locations
file_names=[]
for s in seed_folders:
	file_names.append(os.path.join(dir,os.path.join(s+'/dat/ca.dat')))

def get_data(file_name):
	return np.loadtxt(file_name).T[3]

def array_to_txtfile(arr,file):
    for line in arr:
        file.write(str(line[0])+" "+str(line[1])+" "+str(line[2])+"\n")

#data_shape=np.loadtxt(os.path.join(os.path.join(rxn_path,seed_folders[0]),file_names[0])).shape
#print data_shape
#start=time.time()
times=np.loadtxt(file_names[0]).T[0]

data=[]
for file_name in file_names:	
	data.append(get_data(file_name))
avg=np.average(data,axis=0)
std=np.std(data,axis=0)
#concentration in millimoles
avg*=1e15/N_avo/vol_er
std*=1e15/N_avo/vol_er
f=open(os.path.join(avg_path,"avg.dat"),'w')
array_to_txtfile(np.stack([times,avg,std],axis=-1),f)
f.close()

#p = Pool(processes=8)
#p.map(evaluate, file_names)

#end=time.time()
#print (end-start)

