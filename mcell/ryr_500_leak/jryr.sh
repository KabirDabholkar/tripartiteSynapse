#!/bin/bash

#PBS -p 1
#PBS -j oe
#PBS -J 1-1000:1


/apps/bin/mcell /home/nishant/ryr/2xRSI20V80d${I}.mdl -seed ${PBS_ARRAY_INDEX}
