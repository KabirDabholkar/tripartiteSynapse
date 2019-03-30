#!/usr/bin/python
#PBS -p 1 
#PBS -j oe
#PBS -J 1-1000

import os

# get seed value from array index
seed = os.getenv('PBS_ARRAY_INDEX')

# define bash command
query = '/apps/bin/mcell '+os.getenv('I')+' -seed ' + seed
print query

#print os.getenv('I')


# run bash command
os.system(query)
