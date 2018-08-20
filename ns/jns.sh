#!/bin/bash

#PBS -p 1
#PBS -j oe
#PBS -J 1-2000:1

/apps/bin/mcell /home/nishant/ns/2xNSI20V50d${I}.mdl -seed ${PBS_ARRAY_INDEX}
