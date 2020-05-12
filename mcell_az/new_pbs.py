#!/usr/bin/python
#PBS -p -1
#PBS -j oe
#PBS -J 0-4

import os

# get seed value from array index
array_ind = os.getenv('PBS_ARRAY_INDEX')

file_names=["/home/subhadra/kabir/tripartiteSynapse/mcell_az/simulations/template.mdl"]*5
seeds=[1,2,3,4,5]
seed = str(seeds[int(array_ind)])
file_name = file_names[int(array_ind)]
# define bash command
query = '/apps/bin/mcell '+file_name+' -seed ' + seed
print query

#print os.getenv('I')


# run bash command
os.system(query)
