#!/usr/bin/python
#PBS -p -1
#PBS -J 1-325500

import os
from itertools import product

# get seed value from array index
array_ind = os.getenv('PBS_ARRAY_INDEX')

##############################
#Figure out sim details

def seed_num(vdcc_num):
    #vdcc_num=int(fname.split('V')[1].replace('.mdl',''))
    if vdcc_num>=110 and vdcc_num<=160:
        return 1000
    elif vdcc_num>=90 and vdcc_num<=100:
        return 2000
    elif vdcc_num>=70 and vdcc_num<=80:
        return 3000
    elif vdcc_num>=40 and vdcc_num<=60:
        return 5000
    else:
        return 10

loc="/home/subhadra/kabir/tripartiteSynapse/mcell_az/simulations/"

isi=range(20,181,40)
vdcc=range(60,201,20)
sims=["R150control","R150ER2x","R300ER2x","R150ER3x","R300ER3x"]#,"stores_blocked"]

total_sims=0
for s,i,v in product(sims,isi,vdcc):
    total_sims+=seed_num(v)
    if total_sims>=int(array_ind):
        break
    
    
file_name = loc + s+"_i%d"%i+"v%d"%v
seed = str(int(array_ind)-(total_sims-seed_num(v)))
##############################################


# define bash command
query = '/apps/bin/mcell '+file_name+' -seed ' + seed
print query

#print os.getenv('I')


# run bash command
os.system(query)