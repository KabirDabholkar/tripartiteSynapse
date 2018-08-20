#!/bin/bash

#PBS -p -1
#PBS -j oe
#PBS -J 1-3000:1

/apps/bin/mcell /home/nishant/irs/IRSI40V${I}.mdl -seed ${PBS_ARRAY_INDEX}
